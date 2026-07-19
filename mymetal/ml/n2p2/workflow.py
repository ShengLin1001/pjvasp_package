"""
mymetal.ml.n2p2.workflow

High-level orchestrator for the n2p2 Behler-Parrinello workflow targeting
multi-phase Au (FCC/HCP) potentials. The ``PeiN2p2`` class chains the full
pipeline on top of the lower-level helpers in this subpackage:

    generate_data -> generate_lsf -> select_sf_by_cur -> submit_train
    -> post_training -> check_interface -> post_properties -> post_epoch_scan
    -> post_training_summary / post_check_interface_summary / post_epoch_scan_summary

The three ``post_*_summary`` steps are the only cross-run ones. They re-read the
per-run tables under ``y_post/<run-id>/`` and reduce the training ensemble to a
mean plus a min-max band (or mean +/- std), writing ``*_summary.txt/pdf`` beside
the run directories in ``y_post/``. They never touch a run's own outputs.

It does not hold any model weights itself; it drives the n2p2 tool-chain
(nnp-scaling / nnp-train / nnp-predict) and the LAMMPS ``pair_style hdnnp``
property tests, delegating computation to ``calculate.{sf,cur,post}``,
``dataset.nnpdata``, ``universal.plot.n2p2`` and ``slurm.submit``.

Machine-specific paths (``dir_lmp_utils``, ``path_python``) have no built-in
defaults and must be supplied by the calling project driver script.

Classes:
    PeiN2p2: End-to-end n2p2 training/testing workflow orchestrator.
"""

from mymetal.ml.n2p2.calculate.sf import mysfparams, get_leta_lshift_from_N, generate_g2_g3_g9_blocks, get_largest_rc_from_input_nn
from mymetal.ml.n2p2.calculate.cur import collect_sf_features, filter_zero_columns, cur_select
from mymetal.ml.n2p2.calculate.post import (read_learning_curve, read_normalization, read_trainpoints,
                                            read_trainforces, rmse_me, build_rmse_by_tag_df)
from mymetal.universal.plot.n2p2 import (my_plot_learning_curve, my_plot_compare, my_plot_rmse_by_tag,
                                         my_plot_epoch_stretch, my_plot_epoch_cij, my_plot_epoch_gsfe,
                                         my_plot_epoch_rmse, my_plot_check_interface, VERDICT_CODE)
from mymetal.ml.n2p2.dataset import nnpdata, read_dft_reference
from mymetal.io.general import general_write
from mymetal.slurm.submit import pei_slurm_univ_submit
import numpy as np
import pandas as pd
import os
from datetime import datetime
from pathlib import Path
import shutil
import shlex
import subprocess
import re
import time
import getpass
import warnings


# epoch 扫描的物性字典与列顺序：post_epoch_scan（逐 run 读 LAMMPS 结果）与
# post_epoch_scan_summary（跨 run 汇总）共用同一份定义，保证两者列名/列序一致，
# 也与 mymetal.universal.plot.n2p2 的面板行列顺序对齐。
_LATS = ['fcc', 'bcc', 'hcp']
_CIJ_KEYS = ['C11', 'C12', 'C13', 'C33', 'C44']
# 滑移系：phase 子目录 -> 该相下的 gsfe 类型（与 data_tag_dict 中 A21-2/A22-2 一致）
# FCC 顺序 111 在前、100 在后，与 my_plot_epoch_gsfe 的列顺序保持一致（读取顺序无关）
_GSFE_TYPES = {'fcc': ['FCC_111', 'FCC_100'],
               'hcp': ['HCP_basal', 'HCP_prism1w', 'HCP_pyr1w', 'HCP_pyr2']}
_ALL_GSFE_TYPES = [t for ltype in _GSFE_TYPES.values() for t in ltype]

# c_fcc/c_bcc 立方相 c==a，但绘图 3x3 的 c 行需要，故一并保留
_COLS_STRETCH = ['epoch', 'a_fcc', 'a_bcc', 'a_hcp', 'c_fcc', 'c_bcc', 'c_hcp', 'ca_hcp',
                 'E_fcc', 'E_bcc', 'E_hcp', 'dE_hcp_fcc', 'dE_bcc_fcc']
_COLS_CIJ = ['epoch'] + [f'{k}_{lat}' for lat in _LATS for k in _CIJ_KEYS]
_COLS_GSFE = ['epoch'] + [f'{p}_{t}' for t in _ALL_GSFE_TYPES for p in ['usf', 'sf']]

# 三类 epoch 扫描产物的共同元数据：kind -> (列顺序, 表格浮点格式, 表头量纲说明)
# post_epoch_scan_summary 据此遍历三类产物，避免为每类各写一遍相同的读表/统计/写表逻辑。
_EPOCH_KINDS = {
    'stretch': (_COLS_STRETCH, '14.6f',
                'FCC/BCC a/c are conventional cubic lengths; HCP a/c are hexagonal; '
                'E (eV/atom), dE (meV/atom)'),
    'cij': (_COLS_CIJ, '12.4f',
            'Cij (GPa); cubic phases use C11=C33, C12=C13'),
    'gsfe': (_COLS_GSFE, '12.4f',
             'usf = max(gamma), sf = last gamma value; units mJ/m^2'),
}

# 训练误差表的列（post_training 写出，post_training_summary 读回）：
#   _COLS_RMSE            p_post_learning_curve.txt：逐 epoch 的训练集 E/F RMSE
#   _COLS_RMSE_BY_TAG     p_post_rmse_by_tag.txt 中需要跨 run 统计的误差列
#   _COLS_RMSE_BY_TAG_COUNT  同表中的计数列（数据集属性，各 run 相同，不做统计、原样透传）
_COLS_RMSE = ['epoch', 'E_RMSE_meV/at', 'F_RMSE_meV/A']
_COLS_RMSE_BY_TAG = ['tag', 'E_RMSE_meV/at', 'E_ME_meV/at', 'F_RMSE_meV/A', 'F_ME_meV/A']
_COLS_RMSE_BY_TAG_COUNT = ['n_struct', 'n_fcomp']

# post_*_summary 的色带/误差棒取法 -> 图例文案（{n} 由 run 数填充）。
# 'minmax' 是逐点的跨 run 上下限；'std' 是 mean ± 样本标准差（chap_4 图 1.3 的透明背景读法）。
_BAND_LABEL = {'minmax': 'Min-max over {n} runs',
               'std': 'Mean $\\pm$ std over {n} runs'}

# post_epoch_scan_summary 写出的缺失清单 p_post_absence.txt：列顺序 + status 取值含义。
# 汇总用 nan-aware 统计把缺失静默吃掉（只在屏幕上告警），这张表把它落盘成可追查的记录。
_COLS_ABSENCE = ['run', 'kind', 'epoch', 'status', 'item']
_ABSENCE_STATUS = {
    'missing_table': 'run 没有该 kind 的 p_post_epoch_<kind>.txt（epoch = -1, item = 表名）',
    'missing_epoch': '表在，但没有这个 epoch 的行 —— post_epoch_scan 当时判定该 epoch 输入不全（item = *）',
    'missing_value': '有该 epoch 的行，但 item 这一列是 NaN/inf',
}


class PeiN2p2:
    """n2p2 势函数训练与物性测试的高层调度器。

    该类负责组织训练数据、对称函数、训练目录、后处理目录和 LAMMPS
    物性测试目录，并调用底层 n2p2、LAMMPS、CUR 和绘图工具完成完整流程。

    Attributes:
        dir_root: 项目根目录，必须为绝对路径。
        dir_data: n2p2 数据目录。
        dir_sf: 对称函数生成、scaling 与 CUR 筛选目录。
        dir_train: 训练输入与训练运行目录。
        dir_file: 通用模板和提交脚本目录。
        dir_lmp_utils: LAMMPS 物性测试模板与后处理工具目录。
        path_python: 物性测试后处理所用 Python 解释器路径。
        lele: 从 input.nn.nosf 读取的元素列表。

    TODO:
        - 将初始训练集的每个atoms的原子数调整到36-72，离48最近
        - 选取最后的势函数作MD升温模拟，挑选出某些结构，放到DFT中计算，并append到训练集，形成迭代式训练
        - 添加收集各个势函数的平衡结构
        - 将平衡结构用DFT再计算，并append到训练集，形成迭代式训练
    """

    def __init__(self, dir_data: Path = Path("./data"), dir_sf: Path = Path("./sf"), dir_train: Path = Path("./train"), dir_file: Path = Path("./file"),
                  dir_root: Path = None, if_clean: bool = False,
                  dir_lmp_utils: Path = None,
                  path_python: Path = None):
        """初始化 n2p2 工作流目录、关键路径和元素信息。

        Args:
            dir_data: 数据目录。
            dir_sf: 对称函数与 scaling 目录。
            dir_train: 训练目录。
            dir_file: 模板文件与提交脚本目录。
            dir_root: 项目根目录，必须为绝对路径。
            if_clean: 是否删除已有的 sf 和 train 目录后重建。
            dir_lmp_utils: LAMMPS 物性测试模板与后处理工具目录。
            path_python: 物性测试后处理使用的 Python 解释器路径。

        Raises:
            ValueError: dir_root 或 path_python 不是绝对路径，或缺少必要机器相关路径。
            FileNotFoundError: 找不到 input.nn.nosf。
        """

        # dir_root 必须是绝对路径，否则后续的路径操作会出问题
        if not Path(dir_root).is_absolute():
            raise ValueError(f"dir_root must be an absolute path, got: {dir_root}")
        else:
            self.dir_root = Path(dir_root)

        self.dir_data = Path(dir_data)
        self.dir_sf = Path(dir_sf)
        self.dir_train = Path(dir_train)
        self.dir_file = Path(dir_file)
        # lmp_utils 唯一可信源（LAMMPS 物性测试的 template/post/sub_lmp.bash），
        # properties 方法从这里整文件夹拷成每个 run 的自包含快照。
        # 本机专属路径不进通用包，无内置默认，必须由调用方（项目驱动脚本）显式传入
        if dir_lmp_utils is None:
            raise ValueError("dir_lmp_utils must be specified (no built-in default); "
                             "pass the lmp_utils source dir explicitly.")
        self.dir_lmp_utils = Path(dir_lmp_utils)

        if path_python is None:
            raise ValueError("path_python must be specified (no built-in default); "
                             "pass the python interpreter path explicitly.")
        if not Path(path_python).is_absolute():
            raise ValueError(f"path_python must be an absolute path, got: {path_python}")
        else:
            self.path_python = Path(path_python)

        self.path_input_nn_nosf = self.dir_file / 'input.nn.nosf'
        self.path_input_nn_allsf = self.dir_sf / 'file' / 'input.nn.allsf'
        self.path_input_nn_selectedsf = self.dir_sf / 'file' / 'input.nn.selectedsf'

        # 除了data之外，其他都可以全部删除。
        # if_clean=False 时不删除，避免重跑脚本时误删 sf/y_scaling 的 nnp-scaling 结果
        if if_clean:
            for d in [self.dir_sf, self.dir_train]:
                if d.exists():
                    shutil.rmtree(d)

        os.makedirs(self.dir_data / 'train', exist_ok=True)
        os.makedirs(self.dir_data / 'test', exist_ok=True)
        os.makedirs(self.dir_sf / 'file', exist_ok=True)
        os.makedirs(self.dir_sf / 'y_scaling' / 'y_dir', exist_ok=True)
        os.makedirs(self.dir_train, exist_ok=True)

        if not self.path_input_nn_nosf.is_file():
            raise FileNotFoundError(f"❌ Missing input.nn.nosf: {self.path_input_nn_nosf}. Please make sure it exists in the train directory.")

        self.lele = self._read_elements_from_input_nn_nosf()

        # 最近一次 sbatch 提交出去的 child 作业号（generate_lsf / submit_train / post_properties
        # 每次提交后刷新）。1 核常驻控制器据此调用 wait_jobs() 阻塞等待该阶段全部 child 完成。
        self.last_jobids: list = []


    def _read_elements_from_input_nn_nosf(self) -> list:
        """从 input.nn.nosf 读取 n2p2 元素列表。

        Returns:
            input.nn.nosf 中声明的元素列表。

        Raises:
            ValueError: 缺少 number_of_elements、elements，或两者数量不一致。
        """
        number_of_elements = None
        elements = None

        with open(self.path_input_nn_nosf, "r", encoding="utf-8") as f:
            for line in f:
                # 把当前行按照第一个 # 分割，只保留 # 前面的内容
                fields = line.split("#", 1)[0].split()
                if not fields:
                    continue

                if fields[0] == "number_of_elements":
                    if len(fields) < 2:
                        raise ValueError(f"❌ Missing value after number_of_elements in {self.path_input_nn_nosf}.")
                    number_of_elements = int(fields[1])
                elif fields[0] == "elements":
                    if len(fields) < 2:
                        raise ValueError(f"❌ Missing values after elements in {self.path_input_nn_nosf}.")
                    elements = fields[1:]

        if number_of_elements is None:
            raise ValueError(f"❌ Missing number_of_elements in {self.path_input_nn_nosf}.")
        if elements is None:
            raise ValueError(f"❌ Missing elements in {self.path_input_nn_nosf}.")
        if len(elements) != number_of_elements:
            raise ValueError(
                f"❌ Element count mismatch in {self.path_input_nn_nosf}: "
                f"number_of_elements={number_of_elements}, elements={elements}."
            )

        return elements


    @staticmethod
    def _select_symmetric_strain_dirs(ldir: list, n_side: int = 10, zero_atol: float = 1e-12) -> list:
        """对应变扫描子目录做对称抽样：保留 0 应变 + 正/负各 n_side 个近似均匀点。

        子目录名形如 ``s_+0.0000`` / ``s_-0.0025``（末尾带符号的应变值）。解析应变后按
        零 / 负 / 正 分三组，各组按应变升序排列，再用 ``np.linspace`` 取含首尾的 n_side
        个近似均匀索引——正负两侧各 n_side 个、中间 1 个零应变，共 ``2*n_side+1`` 个。
        某侧点数不足 n_side 时取该侧全部（去重、不补齐）。

        Note:
            "均匀"指在该侧应变的**索引**上均匀（linspace + round），当原网格等间距时
            即应变均匀；含首尾点，故会保留最大幅度应变（如 ±0.12）。原网格非等距时退化
            为按序号均匀。无法从目录名解析出应变的子目录会被跳过。

        Args:
            ldir: 应变扫描的结构子目录列表（每个是一个应变点）。
            n_side: 正负两侧各保留的点数。
            zero_atol: 判定"零应变"的容差。

        Returns:
            抽样后的子目录列表，按应变升序（负 -> 0 -> 正）。
        """
        def _parse_strain(name):
            # 末尾带符号的浮点即应变；解析不到（非应变目录）返回 None，交由调用处跳过
            m = re.search(r'([-+]?\d+(?:\.\d+)?)\s*$', name)
            return float(m.group(1)) if m else None

        def _even_pick(lseq, n):
            # 索引上取 n 个含首尾的近似均匀点；点数不足则全取（round 后去重防重复索引）
            if len(lseq) <= n:
                return list(lseq)
            lidx = np.unique(np.linspace(0, len(lseq) - 1, n).round().astype(int))
            return [lseq[i] for i in lidx]

        lstrain = [(_parse_strain(d.name), d) for d in ldir]
        lstrain = sorted(((s, d) for s, d in lstrain if s is not None), key=lambda x: x[0])
        lzero = [d for s, d in lstrain if abs(s) <= zero_atol]
        lneg  = [d for s, d in lstrain if s < -zero_atol]   # 升序：最负 -> 最接近 0
        lpos  = [d for s, d in lstrain if s >  zero_atol]    # 升序：最接近 0 -> 最正

        return _even_pick(lneg, n_side) + lzero + _even_pick(lpos, n_side)


    def generate_data(self, dir_data_source: Path = None, data_tag_dict: dict = None,
                      if_adjust_size: bool = False, size_top: int = 72,
                      size_bottom: int = 36, size_close: int = 48,
                      dict_subsample: dict = None):
        """从 VASP OUTCAR 生成 n2p2 训练数据。

        遍历 ``dir_data_source/<tag>/<subdir>/y_dir`` 下的每个子目录，用
        ``nnpdata.load_from_outcar`` 读末帧，追加写入 ``data/train/`` 下的
        ``input.data``（总集）和 ``input_<tag>.data``（按 tag 分文件），同时把所有
        力 / 能量分别落盘为 ``forces.npy`` / ``energy.npy`` 备查。

        Note:
            当前版本只往训练集（train/）写，测试集（test/）暂留作占位、未使用；
            写入是 append 语义，但每次运行会先清空 train/ 下旧的 input*.data 再重写，
            故重跑安全、不会重复累加。

        Args:
            dir_data_source: 按 tag 分类存放 VASP 计算结果的根目录。
            data_tag_dict: tag 到子目录列表的映射。
            if_adjust_size: 是否把每个结构的原子数规整到 [size_bottom, size_top] 区间
                内并尽量逼近 size_close（仅面内 a,b 超胞复制，nz 恒为 1）。默认 False，
                保持原结构不变。
            size_top: 原子数上限。原子数已超过此值的结构（如大平板）保持原样不复制。
            size_bottom: 原子数下限。
            size_close: 规整目标原子数，复制方案在区间内尽量逼近此值。
            dict_subsample: 可选的按子目录抽样表 ``{子目录路径子串: n_side}``。某个
                ``subdir`` 命中任一 key（子串匹配）时，只加载该子目录下 0 应变 + 正/负各
                ``n_side`` 个近似均匀应变点（见 ``_select_symmetric_strain_dirs``），用于给
                应变点过密的数据（如 HOEC 每模式 97 点）减重；默认 None 表示全部加载。

        Raises:
            FileNotFoundError: 某个数据子目录下缺少 y_dir。
        """
        os.chdir(self.dir_root)
        workdir = self.dir_data
        path_data_train = workdir / 'train'
        path_data_test  = workdir / 'test' # 暂时没用到
        dir_data_source = Path(dir_data_source)

        # Important
        # Clean old dataset files
        for file in path_data_train.iterdir():
            if file.is_file() and file.name.startswith("input") and file.suffix == ".data":
                file.unlink()

        lforces = []
        lenergies = []
        for tag, lsubdir in data_tag_dict.items():
            for subdir in lsubdir:
                path_subdir = Path(os.path.join(dir_data_source, tag, subdir))
                path_ydir   = path_subdir / 'y_dir'

                if not path_ydir.is_dir():
                    raise FileNotFoundError(f"❌ Missing y_dir: {path_ydir}")

                print(f'========{tag} - {subdir}')
                # search in sub-subdir
                ldir_struct = [d for d in sorted(path_ydir.iterdir()) if d.is_dir()]

                # 可选：对应变点过密的子目录（如 HOEC）对称抽样，0 应变 + 每侧 n_side 个均匀点
                n_side = next((n for token, n in (dict_subsample or {}).items() if token in subdir), None)
                if n_side is not None:
                    n_before = len(ldir_struct)
                    ldir_struct = self._select_symmetric_strain_dirs(ldir_struct, n_side=n_side)
                    print(f"  ↘ subsample {tag}/{subdir}: {n_before} -> {len(ldir_struct)} "
                          f"structures (0 应变 + 每侧 {n_side} 个均匀应变点)")

                for d in ldir_struct:
                    if d.is_dir():
                        comment_file = '/'.join(d.parts[8:])
                        mynnpdata = nnpdata()
                        mynnpdata.load_from_outcar(outcarfile=d / 'OUTCAR', index='-1', tag=tag, comment_file=comment_file)
                        if if_adjust_size:
                            mynnpdata.adjust_size(size_top=size_top, size_bottom=size_bottom, size_close=size_close)
                        mynnpdata.write(outfile_name=path_data_train / f'input_{tag}.data', append=True)
                        mynnpdata.write(outfile_name=path_data_train / 'input.data', append=True)
                        mydict = mynnpdata.get_dict()
                        # mydict['lforces']: 150 * 120 * 3, 150 structures
                        # after vstack, (150*120) * 3
                        # many frames
                        lforces.append(np.vstack(mydict['lforces']))
                        lenergies += mydict['lenergies']

        lforces = np.vstack(lforces)
        lenergies = np.array(lenergies)

        np.save(path_data_train / 'forces.npy', lforces)
        np.save(path_data_train / 'energy.npy', lenergies)


    def generate_lsf(self, lrc_dict: dict = None, n_dict: dict = None, lrs_dict: dict = None, llambd: list = None, lzeta: list = None,
                     if_save: bool = False, if_sbatch: bool = False, lfile: list = ['input.data', 'sub.n2p2.scaling.sh']):
        """生成候选对称函数并准备逐函数 nnp-scaling 目录。

        用 ``generate_g2_g3_g9_blocks`` 生成 G2/G3/G9 候选对称函数，写出
        ``input.nn.allsf``（= input.nn.nosf + 全部候选 SF）和 ``SFs_all.dat``；再把
        每条 SF 单独追加到 ``sf/y_scaling/y_dir/<idx>/input.nn`` 后分发，每个子目录额外
        补一条 rc=0.1 的占位 SF（否则只含 1 条 SF 时 n2p2 会报错），供后续逐函数
        ``nnp-scaling`` 收集 function.data。

        Args:
            lrc_dict: 各元素或元素组合的截断半径设置。
            n_dict: 各类对称函数的数量设置。
            lrs_dict: 径向位移参数设置。
            llambd: 角向函数的 lambda 参数列表。
            lzeta: 角向函数的 zeta 参数列表。
            if_save: 是否写出候选对称函数文件并创建 scaling 子目录。
            if_sbatch: 是否提交 nnp-scaling Slurm 作业。
            lfile: 复制到 scaling 目录的输入数据和提交脚本文件名。

        Returns:
            生成的对称函数文本块列表。
        """

        os.chdir(self.dir_root)

        lblock = generate_g2_g3_g9_blocks(lrc_dict=lrc_dict, n_dict=n_dict, lrs_dict=lrs_dict, llambd=llambd, lzeta=lzeta, if_print=False)

        if if_save:
            path_SF_all_include_comments = self.dir_sf / 'file' / 'SFs_all_include_comments.dat'
            path_SF_all = self.dir_sf / 'file' / 'SFs_all.dat'
            path_input_nn_nosf  = self.dir_sf / 'file' / 'input.nn.nosf'
            path_input_nn_allsf = self.dir_sf / 'file' / 'input.nn.allsf'

            shutil.copy(self.path_input_nn_nosf, path_input_nn_nosf)
            shutil.copy(self.path_input_nn_nosf, path_input_nn_allsf)

            with open(path_SF_all_include_comments, "w", encoding="utf-8") as f:
                for block in lblock:
                    f.write(block)

            with open(path_input_nn_allsf, "a", encoding="utf-8") as f:
                for block in lblock:
                    f.write(block)

            with open(path_SF_all, "w", encoding="utf-8") as f:
                num_sf = 0
                for block in lblock:
                    for b in block.strip('\n').split('\n'):
                        # 排除空行和注释行
                        if b.strip() and not b.lstrip().startswith("#"):
                            f.write(b + '\n')
                            num_sf += 1
                print(f"Has {num_sf} sfs.")

            i = 0

            with open(path_SF_all, 'r') as f:
                lines = f.readlines()

            for line in lines:
                if line.strip() and not line.lstrip().startswith("#"):
                    path_subdir = self.dir_sf / 'y_scaling' / 'y_dir' / f'{i:04d}'
                    os.makedirs(path_subdir, exist_ok=True)
                    i += 1

                    shutil.copy(path_input_nn_nosf, path_subdir / 'input.nn')
                    with open(path_subdir / 'input.nn','a') as f:
                        f.write(line)

                        # append two 0 value SF (rc = 0.1), ortherwise n2p2 will break
                        f.write(f'symfunction_short {self.lele[0]} 3 {self.lele[0]} {self.lele[0]} 1  1.000 1.000 0.100\n')

                    # Copy input.data, sub.*
                    shutil.copy(self.dir_data / 'train' / lfile[0], path_subdir / lfile[0])

                    #os.chdir(path_subdir)
                    #print(os.getcwd())
                    #os.system('sbatch sub.*')
                    #os.chdir(self.dir_root)

            shutil.copy(self.dir_file / lfile[1], self.dir_sf / 'y_scaling' / lfile[1])
            self.last_jobids = []
            if if_sbatch:
                # 捕获 nnp-scaling 作业号，供控制器阻塞等待其完成后再做 CUR 选 SF
                self.last_jobids = [self._sbatch_capture(lfile[1], cwd=self.dir_sf / 'y_scaling')]

        return lblock


    def select_sf_by_cur(self, n_select: int = 48, max_zero_frac: float = 1.0, zero_atol: float = 0.0,
                         min_std: float = 1e-4, standardize: bool = True,
                         if_copy_to_train: bool = True, n_select_by_type: dict = None) -> list:
        """用零值/方差过滤和 CUR 从候选对称函数中选择最终集合。

        从每个 scaling 子目录的 function.data 收集每条候选 SF 的逐原子特征，先做一道宽松
        过滤——默认**关掉零值过滤**（``max_zero_frac=1.0``），只用标准差 <= ``min_std``
        （默认 1e-4）剔除近常数（对数据集不敏感）的死特征——再用 CUR 从剩余 SF 中选出最重要
        的 ``n_select`` 条，追加到干净的 input.nn.nosf 后生成最终 input.nn——不含占位 SF，
        因为占位 SF 只写进各 scaling 子目录、从不写进 file/ 的 nosf。

        CUR 直接在逐原子特征矩阵 ``feat_atom``（本例约 8.5 万原子 × 220 SF，单次 ~80s）上
        做，保留帧内原子间的区分度。分工是：过滤只负责剔除真正退化（近常数）的列，避免
        标准化把它们的数值噪声放大到与真信号同量级；真正"选出有用 SF"由 CUR + 列标准化完成。

        Note:
            零值比例只衡量"这条 SF 对多少原子为 0"（稀疏性），而标准差衡量"它在数据集上变不变"
            （敏感度）。"对数据集不敏感"的本意是**不变**，故以 ``min_std`` 为主判据；零值过滤
            默认关闭，因为在本例中它会误杀稀疏但有区分力的角度型 SF（如仅 5.2% 为 0、但满量程
            的 G2/G9）。需要时可通过 ``max_zero_frac<1`` 重新开启。
            ``zero_atol=0.0`` 与参考模板一致（只算精确零）；取 ~1e-9 会把截断球边缘的
            "准零"也算上（n2p2 的 function.data 为 10 位小数）。

        Args:
            n_select: CUR 选择的对称函数数量（仅在 n_select_by_type 为 None 时生效）。
            max_zero_frac: 允许的最大零值比例。默认 1.0 = 关闭零值过滤；设 <1 可重新开启。
            zero_atol: 判定近零值的绝对容差（仅在零值过滤开启时有意义）。
            min_std: 保留 SF 所需的最小逐原子标准差，剔除近常数（对数据集不敏感）的 SF。
                默认 1e-4，作为主过滤判据；这是"误杀稀疏有用角度型"与"放进噪声放大的近平
                SF"之间的折中，需要更激进/更保守可在 1e-3 ~ 1e-6 间调。
            standardize: CUR 前是否对各列做 z-score 标准化，去除未中心化 CUR 偏爱高幅值
                常数列的偏差，使选择跟随变化量而非绝对水平。
            if_copy_to_train: 是否将最终 input.nn 复制到训练目录。
            n_select_by_type: 若给定（{symfunction_type: 数量}，如 {2: 15, 9: 5}），对每个
                对称函数类型单独 CUR、各选够指定数量，忽略 n_select。用于保证类型配额
                （纯体积拉伸下角度型 G9 方差低，全局 CUR 会把它们全部忽略）。若某类型过滤后
                的候选数不足配额，则取二者较小值并打印告警，而不是报错——数量不足恰说明该
                数据集对此类型不敏感。

        Returns:
            被选中的对称函数行列表。

        Raises:
            ValueError: 候选对称函数数量与 scaling 子目录数量不一致。
            FileNotFoundError: 某个 scaling 子目录缺少 function.data。
        """

        os.chdir(self.dir_root)

        dir_sf_file = self.dir_sf / 'file'
        path_SF_all = dir_sf_file / 'SFs_all.dat'
        dir_ydir = self.dir_sf / 'y_scaling' / 'y_dir'

        with open(path_SF_all, 'r', encoding='utf-8') as f:
            lsf_all = [line for line in f if line.strip()]

        # 校验 scaling 已全部完成
        lsubdir = sorted([d for d in dir_ydir.iterdir() if d.is_dir()])
        if len(lsubdir) != len(lsf_all):
            raise ValueError(f"❌ SF count mismatch: {len(lsf_all)} sfs in {path_SF_all}, but {len(lsubdir)} dirs in {dir_ydir}.")
        for d in lsubdir:
            if not (d / 'function.data').is_file():
                raise FileNotFoundError(f"❌ Missing function.data in {d}. Please make sure nnp-scaling finished.")

        # 收集特征：feat_atom (总原子数, n_sf)，feat_av (帧数, n_sf)
        feat_atom, feat_av = collect_sf_features(lsubdir)
        np.save(dir_sf_file / 'feat_atom.npy', feat_atom)
        np.save(dir_sf_file / 'feat_av.npy', feat_av)

        # 剔除零值比例过高或近常数（对数据集不敏感）的 sf
        kept_idx, dropped_idx = filter_zero_columns(feat_atom, max_zero_frac=max_zero_frac,
                                                    zero_atol=zero_atol, min_std=min_std)
        with open(dir_sf_file / 'SFs_dropped.dat', 'w', encoding='utf-8') as f:
            for i in dropped_idx:
                f.write(lsf_all[i])
        print(f"Dropped {len(dropped_idx)} sfs (zero fraction > {max_zero_frac} or std <= {min_std}), "
              f"{len(kept_idx)} sfs left for CUR.")

        # CUR 选择（基于逐原子特征 feat_atom + 列标准化）
        if n_select_by_type is None:
            lidx_selected = cur_select(feat_atom[:, kept_idx], n_select=n_select, standardize=standardize)
            # 返回的lidx_selected是去0后的索引，需要映射回原始索引
            # sort 仅为了美观
            lsf_selected = sorted([lsf_all[kept_idx[i]] for i in lidx_selected])
        else:
            # 分类型配额：对每个 symfunction type 单独 CUR，各选够 n_select_by_type[type]。
            # type 取每行第 3 个字段（symfunction_short <center> <type> ...）。
            def _sf_type(line):
                fields = line.split()
                return int(fields[2]) if len(fields) >= 3 and fields[0] == 'symfunction_short' else None
            lsf_selected = []
            for sftype, n_sel_t in n_select_by_type.items():
                kept_t = [gi for gi in kept_idx if _sf_type(lsf_all[gi]) == int(sftype)]
                # 过滤后候选不足配额：取较小值并告警，而不是报错。数量不足恰恰说明
                # 该数据集对此类型不敏感，强行凑数只会引入近常数死特征。
                n_take = min(n_sel_t, len(kept_t))
                if n_take == 0:
                    print(f"⚠️  type {sftype}: 0 candidate sfs left after zero/std filter, skipped "
                          f"(数据集对该类型完全不敏感).")
                    continue
                if n_take < n_sel_t:
                    print(f"⚠️  type {sftype}: only {len(kept_t)} candidate sfs left after zero/std "
                          f"filter, selecting {n_take}/{n_sel_t} (数据集对该类型敏感的 SF 不足).")
                sel_local = cur_select(feat_atom[:, kept_t], n_select=n_take, standardize=standardize)
                lsf_selected += [lsf_all[kept_t[i]] for i in sel_local]
            lsf_selected = sorted(lsf_selected)
        with open(dir_sf_file / 'SFs_selected.dat', 'w', encoding='utf-8') as f:
            for line in lsf_selected:
                f.write(line)
        by_type_msg = '' if n_select_by_type is None else f' (by type {dict(n_select_by_type)})'
        print(f"Selected {len(lsf_selected)} sfs by CUR{by_type_msg}.")

        # 最终 input.nn = 干净的 input.nn.nosf（不含 dummy sf）+ 选中的 sf
        # 为了保证不含dummy sf, 已经把前面的逻辑改为在subdir里追加sf，
        # 而不是在file里追加，这样就不会把dummy sf也加到input.nn里了
        path_input_nn_selectedsf = dir_sf_file / 'input.nn.selectedsf'
        shutil.copy(self.path_input_nn_nosf, path_input_nn_selectedsf)
        with open(path_input_nn_selectedsf, 'a', encoding='utf-8') as f:
            for line in lsf_selected:
                f.write(line)

        if if_copy_to_train:
            shutil.copy(path_input_nn_selectedsf, self.dir_train / 'input.nn')

        return lsf_selected


    # ===== 1 核常驻控制器的 child 作业编排（提交即捕获作业号，阶段末统一阻塞等待）=====
    @staticmethod
    def _sbatch_capture(script: str, cwd: Path, retries: int = 99, retry_interval: int = 10) -> str:
        """在 cwd 下 sbatch 单个脚本并返回作业号（--parsable，确定性捕获）。

        提交偶发失败（slurmctld 繁忙 / 超时 / 通信抖动）时 returncode 非 0、拿不到作业号。
        这类失败意味着作业根本没被创建，重试绝对安全（不会产生重复作业），故轮询重试而不是
        立刻抛错中断整条流水线；重试用尽仍失败才抛 RuntimeError。与 _snapshot_jobids 的
        squeue 重试同一套 99×10s 策略。

        Args:
            script: 提交脚本文件名（相对 cwd）。
            cwd: 提交目录（作业的工作目录即此目录）。
            retries: 提交失败时的重试次数。
            retry_interval: 每次重试前等待的秒数。

        Returns:
            Slurm 作业号字符串。

        Raises:
            RuntimeError: 重试用尽仍无法提交，或提交成功却未解析到作业号。
        """
        last_err = ''
        for attempt in range(1, retries + 1):
            out = subprocess.run(['sbatch', '--parsable', str(script)], cwd=str(cwd),
                                 capture_output=True, text=True)
            # returncode==0 才算提交成功：--parsable 无 --wait，提交成功即 0、stdout 是作业号。
            if out.returncode == 0:
                # --parsable 输出 "<jobid>" 或 "<jobid>;<cluster>"
                jobid = out.stdout.strip().split(';')[0]
                if not jobid:
                    raise RuntimeError(f"❌ Could not parse jobid from sbatch in {cwd}: {out.stdout!r}")
                print(f"  submitted job {jobid}: {cwd}/{script}")
                return jobid
            last_err = out.stderr.strip() or out.stdout.strip()
            if attempt < retries:
                print(f"  ⚠️  sbatch submit failed (attempt {attempt}/{retries}) in {cwd}: {last_err}; retry in {retry_interval}s")
                time.sleep(retry_interval)
        raise RuntimeError(f"❌ sbatch failed after {retries} attempts in {cwd}: {last_err}")


    @staticmethod
    def _snapshot_jobids(retries: int = 99, retry_interval: int = 10) -> set:
        """当前用户排队/运行中的全部作业号集合。

        用于捕获经第三方提交器（如 pei_slurm_univ_submit）提交、无法直接拿到作业号的
        child 作业：提交前后各取一次快照，差集即本次新提交的作业。

        关键：squeue 偶发失败（slurmctld 繁忙/超时）时，stdout 为空。**绝不能**把空输出
        当作"队列里没作业"——那会让 wait_jobs 误判作业已离队、在训练还在跑时就开始后处理
        （曾导致读到半截 learning-curve 的非 5 倍数 epoch、缺权重文件而崩溃）。故这里检查
        returncode，失败则重试，重试用尽仍失败就抛 RuntimeError 让调用方显式处理，而不是
        静默返回空集。

        Args:
            retries: squeue 失败时的重试次数。
            retry_interval: 每次重试前等待的秒数。

        Returns:
            作业号字符串集合。

        Raises:
            RuntimeError: squeue 连续失败、无法获得可靠的队列快照。
        """
        user = os.environ.get('USER') or getpass.getuser()
        last_err = ''
        for attempt in range(1, retries + 1):
            try:
                out = subprocess.run(['squeue', '-u', user, '-h', '-o', '%i'],
                                     capture_output=True, text=True)
            except Exception as e:
                # try 发生异常走这里
                last_err = repr(e)
            else:
                # try 成功走这里
                # returncode==0 才信任 stdout；非 0 时空 stdout 不代表队列为空
                if out.returncode == 0:
                    return set(x.strip() for x in out.stdout.split() if x.strip())
                last_err = f'returncode={out.returncode}, stderr={out.stderr.strip()}'
            if attempt < retries:
                print(f"  ⚠️  squeue failed (attempt {attempt}/{retries}): {last_err}; retry in {retry_interval}s")
                time.sleep(retry_interval)
        raise RuntimeError(f"❌ squeue failed after {retries} attempts: {last_err}")


    @staticmethod
    def wait_jobs(ljobid: list, poll_interval: int = 30, label: str = '') -> None:
        """阻塞直到给定作业号全部离开队列（完成/失败/取消均算离队）。

        仅判断是否仍在队列，不区分成功失败；失败由下游产物检查（post_training 等的
        FileNotFoundError）兜底暴露。

        Args:
            ljobid: 需要等待的作业号列表。
            poll_interval: 轮询 squeue 的间隔（秒）。
            label: 打印用的阶段标签。
        """
        ljobid = set(str(j) for j in ljobid if str(j).strip())
        if not ljobid:
            print(f"  [wait:{label}] no jobs to wait for.")
            return
        print(f"  [wait:{label}] waiting for {len(ljobid)} job(s): {sorted(ljobid)}")
        t0 = datetime.now()
        while True:
            try:
                remaining = PeiN2p2._snapshot_jobids() & ljobid
            except RuntimeError as e:
                # 快照不可靠时只能保守地认为作业仍在队列，继续等——
                # 宁可多等，也绝不在拿不到可靠队列状态时误判作业已结束
                print(f"  [wait:{label}] {e}; assume jobs still running, keep waiting.")
                time.sleep(poll_interval)
                continue
            if not remaining:
                break
            time.sleep(poll_interval)
        dt = (datetime.now() - t0).total_seconds()
        print(f"  [wait:{label}] all {len(ljobid)} job(s) left the queue after {dt:.0f}s.")


    def submit_train(self, dir_run: Path = Path('./train/y_n2p2_train/y_dir/001'),
                     lfile: list = ['input.data', 'sub.n2p2.train.sh'], if_sbatch: bool = False,
                     random_seed: int = None, sed_overrides: dict = None):
        """准备单次 n2p2 训练目录并可选提交作业。

        在 ``dir_run`` 下汇齐最终 input.nn（CUR 选出的 SF）、input.data 和提交脚本后
        sbatch。作业内会先用该 input.nn 重新跑 ``nnp-scaling`` 生成 scaling.data
        （nnp-train 必读；逐函数 scaling 阶段每条 SF 单独的 scaling.data 不能复用），
        再跑 ``nnp-train``。为防覆盖已有训练，目标目录若已存在 scaling.data /
        learning-curve.out 会直接报错——重训请换一个 dir_run（如 y_dir/002）。

        Args:
            dir_run: 本次训练运行目录。
            lfile: 训练数据文件和提交脚本文件名。
            if_sbatch: 是否提交 Slurm 作业。
            random_seed: 若给定，则在拷贝 input.nn 后、sbatch 前 sed 改写本 run input.nn
                的 random_seed。多 seed ensemble 由此外部控制（各 dir_run 用不同种子）。
            sed_overrides: 若给定（{关键词: 取值}），在拷贝 input.nn 后、sbatch 前逐个 sed 就地
                改写本 run input.nn 的对应标量关键词（如 {'force_weight': 0.6,
                'short_force_fraction': 0.04}）。参数扫描由此外部控制（各 dir_run 用不同超参、
                同一 random_seed）。改后立即校验，未命中即报错以防扫描退化为同一参数。

        Raises:
            FileNotFoundError: 缺少 input.nn、input.data 或提交脚本。
            FileExistsError: 目标训练目录中已存在训练产物。
            RuntimeError: sed 改写 random_seed / sed_overrides 失败或未生效。
        """

        os.chdir(self.dir_root)
        dir_run = Path(dir_run)

        path_input_nn = self.dir_train / 'input.nn'
        path_input_data = self.dir_data / 'train' / lfile[0]
        path_sub = self.dir_file / lfile[1]

        if not path_input_nn.is_file():
            raise FileNotFoundError(f"❌ Missing {path_input_nn}. Please run select_sf_by_cur first.")
        if not path_input_data.is_file():
            raise FileNotFoundError(f"❌ Missing {path_input_data}. Please run generate_data first.")
        if not path_sub.is_file():
            raise FileNotFoundError(f"❌ Missing submit script: {path_sub}.")

        # 防止覆盖已完成/进行中的训练，重训请换一个 dir_run（如 y_dir/002）
        for fname in ['scaling.data', 'learning-curve.out']:
            if (dir_run / fname).is_file():
                raise FileExistsError(f"❌ {dir_run / fname} already exists. "
                                      f"Use a new dir_run to avoid overwriting a previous training run.")

        os.makedirs(dir_run, exist_ok=True)
        shutil.copy(path_input_nn, dir_run / 'input.nn')
        shutil.copy(path_input_data, dir_run / lfile[0])
        shutil.copy(path_sub, dir_run / lfile[1])

        # 多 seed ensemble：拷贝之后、sbatch 之前用 sed 改写本 run input.nn 的 random_seed，
        # 只换数字、保留对齐与行尾注释（^(random_seed<空白>)<数字>）。
        if random_seed is not None:
            path_run_nn = dir_run / 'input.nn'
            seed = int(random_seed)
            rc = os.system("sed -i -E 's/^(random_seed[[:space:]]+)[0-9]+/\\1%d/' %s"
                           % (seed, shlex.quote(str(path_run_nn))))
            if os.waitstatus_to_exitcode(rc) != 0:
                raise RuntimeError(f"❌ sed random_seed failed for {path_run_nn} (exit {rc}).")
            # 校验：模板缺 random_seed 行时 sed 会静默不改，导致 ensemble 退化为同一 seed
            import re as _re
            m = _re.search(r'(?m)^random_seed\s+(\d+)', path_run_nn.read_text(encoding='utf-8'))
            if m is None or int(m.group(1)) != seed:
                raise RuntimeError(f"❌ random_seed not set to {seed} in {path_run_nn} "
                                   f"(found: {None if m is None else m.group(1)}).")
            print(f"  set random_seed={seed} in {path_run_nn}")

        # 参数扫描：拷贝之后、sbatch 之前用 sed 就地改写本 run input.nn 的任意标量关键词
        # （只换第一个取值 token、保留对齐与行尾注释），改后逐个校验。
        if sed_overrides:
            import re as _re
            path_run_nn = dir_run / 'input.nn'
            for key, val in sed_overrides.items():
                val_str = f'{val:g}' if isinstance(val, float) else str(val)
                rc = os.system("sed -i -E 's/^(%s[[:space:]]+)[^[:space:]#]+/\\1%s/' %s"
                               % (key, val_str, shlex.quote(str(path_run_nn))))
                if os.waitstatus_to_exitcode(rc) != 0:
                    raise RuntimeError(f"❌ sed {key} failed for {path_run_nn} (exit {rc}).")
                m = _re.search(rf'(?m)^{_re.escape(key)}\s+(\S+)',
                               path_run_nn.read_text(encoding='utf-8'))
                if m is None or abs(float(m.group(1)) - float(val)) > 1e-12:
                    raise RuntimeError(f"❌ {key} not set to {val_str} in {path_run_nn} "
                                       f"(found: {None if m is None else m.group(1)}).")
                print(f"  set {key}={val_str} in {path_run_nn}")

        n_sf = 0
        with open(path_input_nn, 'r', encoding='utf-8') as f:
            for line in f:
                fields = line.split('#', 1)[0].split()
                if fields and fields[0] == 'symfunction_short':
                    n_sf += 1
        print(f"Prepared training run in {dir_run} with {n_sf} sfs.")

        self.last_jobids = []
        if if_sbatch:
            # 捕获本 run 的 nnp-train 作业号；控制器把多 seed 的作业号汇总后统一阻塞等待
            self.last_jobids = [self._sbatch_capture(lfile[1], cwd=dir_run)]


    @staticmethod
    def _write_table(path: Path, header: list, df: pd.DataFrame, float_format: str = '16.8f') -> None:
        """写出带注释头的格式化表格。

        先写入 header 注释行，再用通用模板 ``general_write`` 追加对齐格式化的 DataFrame
        （含列名），末尾补一个空行。本类所有 p_post_*.txt 表格都经此输出，以保证统一的
        注释头 + 列对齐风格。

        Args:
            path: 输出文件路径。
            header: 写在表格前的注释行列表。
            df: 待写出的 DataFrame。
            float_format: 浮点数格式。
        """

        with open(path, 'w', encoding='utf-8') as f:
            for line in header:
                f.write(line + '\n')
        general_write(path, if_append=True, dfc=df, if_write_col_num=True, float_format=float_format)
        with open(path, 'a', encoding='utf-8') as f:
            f.write('\n')


    # y_n2p2_train/y_dir
    def post_training(self, dir_run: Path = Path('./train/y_n2p2_train/y_dir/001'), epoch: int = None) -> Path:
        """后处理 n2p2 训练结果并生成误差表格和图像。

        ``y_post/<run-id>/`` 与 ``y_dir/<run-id>`` 的训练一一对应；本方法生成
        ``training/`` 下的全部产物（均为物理单位）：

        - ``p_post_learning_curve.txt/pdf``：epoch vs 能量/力 RMSE（meV/atom、meV/Å）。
        - ``epoch.txt``：选定 epoch（默认末 epoch；可手改后用 ``epoch=`` 重跑）。
        - ``trainpoints/trainforces.data``：选定 epoch 的 DFT-vs-NNP 对照，附 tag。
        - ``p_post_rmse.txt/pdf``：总体 E/F 的 RMSE、ME（表）+ 按 tag 着色的散点图（原 dncompare）。
        - ``p_post_rmse_by_tag.txt/pdf``：按 tag 拆分的能量/力误差（表 + 柱状图）。

        Note:
            训练开了 ``normalize_data_set``，故 trainpoints/trainforces 是训练单位，须用
            nnp-train 写回 input.nn 的 mean_energy/conv_energy/conv_length 换算回物理单位；
            而 learning-curve.out 的 pu 列本就是物理单位（RMSE 是差值统计量，mean_energy
            自然消去）。

        Args:
            dir_run: 训练运行目录。
            epoch: 需要分析的 epoch；为 None 时使用最后一个 epoch。

        Returns:
            本次训练对应的后处理目录。

        Raises:
            FileNotFoundError: 缺少 learning-curve、指定 epoch 的预测文件或权重文件。
            ValueError: input.data 中结构数与 trainpoints 中结构数不一致。
        """
        os.chdir(self.dir_root)
        dir_run = Path(dir_run)

        path_lc = dir_run / 'learning-curve.out'
        if not path_lc.is_file():
            raise FileNotFoundError(f"❌ Missing {path_lc}. Training not finished or wrong dir_run?")

        # 目录骨架：training/ 是本步产物；properties/ 由 post_properties 整体组装，这里不再预建
        dir_post = dir_run.parent.parent / 'y_post' / dir_run.name
        dir_training = dir_post / 'training'
        os.makedirs(dir_training, exist_ok=True)

        # 1. 训练曲线（col2/col10 为物理单位 pu (eV/atom, eV/Å)）
        lc = read_learning_curve(path_lc)
        epochs = lc[:, 0].astype(int)
        e_rmse_mev = lc[:, 1] * 1e3
        f_rmse_mev = lc[:, 9] * 1e3
        df_lc = pd.DataFrame({'epoch': epochs, 'E_RMSE_meV/at': e_rmse_mev, 'F_RMSE_meV/A': f_rmse_mev})
        self._write_table(dir_training / 'p_post_learning_curve.txt',
                          ['# Training-set RMSE per epoch, physical units'], df_lc, float_format='14.6f')
        my_plot_learning_curve(epochs, e_rmse_mev, f_rmse_mev, dir_training / 'p_post_learning_curve.pdf')

        # 2. 选 epoch（默认末 epoch；写出文件的间隔由 input.nn 的 write_weights_epoch 决定）
        if epoch is None:
            epoch = int(epochs[-1])
        path_tp = dir_run / f'trainpoints.{epoch:06d}.out'
        path_tf = dir_run / f'trainforces.{epoch:06d}.out'
        lweights = sorted(dir_run.glob(f'weights.*.{epoch:06d}.out'))
        if not path_tp.is_file() or not path_tf.is_file() or not lweights:
            raise FileNotFoundError(f"❌ Missing trainpoints/trainforces/weights for epoch {epoch} in {dir_run}.")
        with open(dir_training / 'epoch.txt', 'w', encoding='utf-8') as f:
            f.write(f'{epoch:06d}\n')

        # 3. 训练单位 -> 物理单位. Training units -> physical units.
        mean_energy, conv_energy, conv_length = read_normalization(dir_run / 'input.nn')
        tp = read_trainpoints(path_tp)
        tf = read_trainforces(path_tf)
        e_ref = tp[:, 1] / conv_energy + mean_energy
        e_nnp = tp[:, 2] / conv_energy + mean_energy
        f_ref = tf[:, 2] * conv_length / conv_energy
        f_nnp = tf[:, 3] * conv_length / conv_energy

        # tag：用 nnpdata 的轻量方法解析 input.data（单一数据来源、~30x 快于 load_from_datafile），
        # trainpoints/trainforces 第 1 列 index = input.data 中 0 起的结构顺序号（已验证）
        ltag = nnpdata.read_tags_from_datafile(dir_run / 'input.data')
        if len(ltag) != tp.shape[0]:
            raise ValueError(f"❌ Structure count mismatch: {len(ltag)} in input.data, "
                             f"{tp.shape[0]} in {path_tp}.")
        tag_e = np.array([ltag[int(i)] for i in tp[:, 0]])
        tag_f = np.array([ltag[int(i)] for i in tf[:, 0]])

        df_tp = pd.DataFrame({'index': tp[:, 0].astype(int), 'tag': tag_e,
                              'E_dft_eV/at': e_ref, 'E_nnp_eV/at': e_nnp})
        self._write_table(dir_training / 'trainpoints.data',
                          [f'# Epoch {epoch:06d}, physical units'], df_tp, float_format='16.8f')
        df_tf = pd.DataFrame({'index_s': tf[:, 0].astype(int), 'index_a': tf[:, 1].astype(int), 'tag': tag_f,
                              'F_dft_eV/A': f_ref, 'F_nnp_eV/A': f_nnp})
        self._write_table(dir_training / 'trainforces.data',
                          [f'# Epoch {epoch:06d}, physical units'], df_tf, float_format='16.8f')

        # 4. 按 tag 拆分误差（ME = mean(NNP - DFT)），并可视化
        df_tag = build_rmse_by_tag_df(e_ref, e_nnp, tag_e, f_ref, f_nnp, tag_f)
        self._write_table(dir_training / 'p_post_rmse_by_tag.txt',
                          [f'# Epoch {epoch:06d}, physical units. E per structure (per atom), F per component.',
                           '# ME = mean(NNP - DFT)'], df_tag, float_format='12.3f')
        my_plot_rmse_by_tag(df_tag, dir_training / 'p_post_rmse_by_tag.pdf')

        # 5. 总体 RMSE/ME（表 p_post_rmse.txt）+ 能量/力散点图（图 p_post_rmse.pdf，原 dncompare）
        e_rmse, e_me = rmse_me(e_ref, e_nnp)
        f_rmse, f_me = rmse_me(f_ref, f_nnp)
        df_rmse = pd.DataFrame({'quantity': ['E_meV/at', 'F_meV/A'],
                                'RMSE': [e_rmse * 1e3, f_rmse * 1e3],
                                'ME': [e_me * 1e3, f_me * 1e3]})
        self._write_table(dir_training / 'p_post_rmse.txt',
                          [f'# Epoch {epoch:06d}, overall DFT-vs-NNP error, physical units.',
                           '# E per structure (per atom), F per component. ME = mean(NNP - DFT)'],
                          df_rmse, float_format='12.3f')
        my_plot_compare(e_ref, e_nnp, tag_e, f_ref, f_nnp, tag_f, dir_training / 'p_post_rmse.pdf',
                     text_e=f'RMSE = {e_rmse * 1e3:.3f} meV/atom\nME = {e_me * 1e3:.3f} meV/atom',
                     text_f=f'RMSE = {f_rmse * 1e3:.3f} meV/$\\mathrm{{\\AA}}$\nME = {f_me * 1e3:.3f} meV/$\\mathrm{{\\AA}}$')

        print(f"Post-processed training in {dir_training}: epoch {epoch:06d}, "
              f"E_RMSE {e_rmse * 1e3:.3f} meV/atom, F_RMSE {f_rmse * 1e3:.3f} meV/A.")
        return dir_post


    def _prepare_epoch_dir(self, dir_epoch: Path, dir_run: Path, lw: list, dir_lmp_utils: Path) -> None:
        """为单个 epoch 铺好 LAMMPS 物性测试目录（仅脚手架，不做 sed）。

        在 ``<epoch>/`` 下准备三样东西，供 runner 在作业里就地 sed / 运行：

        - ``potential/``：势函数三件套——input.nn + scaling.data（取自训练 run）+
          每元素一份 ``weights.<elem>.data``（hdnnp 按 cwd 的 dir 读，命名同 check_interface）。
        - ``template/``：整目录拷自 ``dir_lmp_utils/template``，runner 据此 sed
          pair_style / pair_coeff 等。
        - ``post/``：整目录拷自 ``dir_lmp_utils/post``，供后处理脚本使用。

        Note:
            历史上 stretch/Cij/gsfe 的逐相 sed（lat、aa/cc、pair_style 等）在 Python 这里做；
            现已全部下推到 runner（pei_lmp_run_properties）在作业里完成，本方法只负责铺目录、
            不再碰模板内容。

        Args:
            dir_epoch: 当前 epoch 的目标目录。
            dir_run: n2p2 训练运行目录。
            lw: 当前 epoch 的权重文件列表。
            dir_lmp_utils: LAMMPS 模板与后处理工具源目录。
        """

        # potential/：势函数三件套（input.nn/scaling.data 取自训练 run；weights 每元素一份，复用 check_interface 同款命名）
        dir_potential = dir_epoch / 'potential'
        dir_template  = dir_epoch / 'template'
        dir_post      = dir_epoch / 'post'
        os.makedirs(dir_potential, exist_ok=True)
        shutil.copy(dir_run / 'input.nn', dir_potential / 'input.nn')
        shutil.copy(dir_run / 'scaling.data', dir_potential / 'scaling.data')
        for w in lw:
            shutil.copy(w, dir_potential / f"weights.{w.name.split('.')[1]}.data")

        shutil.copytree(dir_lmp_utils / "template", dir_template, dirs_exist_ok=True)  # 整目录拷一份到 epoch 里，后续 sed 修改模板里的 pair_style/pair_coeff
        shutil.copytree(dir_lmp_utils / "post", dir_post, dirs_exist_ok=True)          # 整目录拷一份到 epoch 里，供后续提交脚本用（如注入 lmp_utils 源路径）


    @staticmethod
    def _dict_to_cli_args(dict_args: dict) -> str:
        """将参数字典转换为 pei_slurm_univ_submit.py 的 CLI 长选项字符串。

        每个 (key, value) 一律渲染成 ``--key <value>``（CLI 全用 -- 长选项，dest 与库形参
        同名）；值统一 ``str`` 化后用 ``shlex.quote`` 保护，关键是 cmd 里的空格 / 双引号 /
        ``$(pwd)`` 不被本层 shell 解析——原样塞进 CLI、写进生成的作业脚本，留到作业里
        （cwd=epoch 目录）才展开 ``$(pwd)``。

        Note:
            布尔值也按 ``--key True`` / ``--key False`` 输出，不做 store_true 式的省略；
            CLI 侧 ``--if_sbatch`` / ``--if_use_my_launcher`` 用
            ``type=parse_bool, nargs="?", const=True`` 接收，故能正确解析这种带值写法。

        Args:
            dict_args: 参数名到参数值的映射。

        Returns:
            经过 shell 转义的 CLI 参数字符串。
        """
        parts = []
        for key, value in dict_args.items():
            parts.append(f'--{key} {shlex.quote(str(value))}')
        return ' '.join(parts)


    # y_n2p2_train/y_post
    def post_properties(self, dir_run: Path = Path('./train/y_n2p2_train/y_dir/001'),
                        if_sbatch: bool = False,
                        dict_args_to_submit: dict = {'preset': 'zcm6-lammps-0',
                                                     'chunks': 5,
                                                     'nodes': 1,
                                                     'ncores': 1,
                                                     'launcher_type': 'srun',
                                                     'if_use_my_launcher': True,},
                        ) -> Path:
        """准备并可选提交跨 epoch 的 LAMMPS hdnnp 物性测试。

        Python 只搭脚手架，真正的 sed / 运行交给单目录 runner ``pei_lmp_run_properties``，
        由通用引擎 ``pei_slurm_univ_submit`` 提交，一步到位：

        1. 逐 epoch 用 :meth:`_prepare_epoch_dir` 铺好
           ``properties/y_epoch_scan/y_dir/<epoch>/``（potential/ + template/ + post/，不做 sed）。
        2. 生成传给 runner 的 4 个参数：pair_style、pair_coeff、python_path、mass_content
           （第 5 个 lmp_template_path 用 runner 默认的 ./template，不在此传）；所有 sed 都由
           runner 在作业里做。
        3. 每个 ``<epoch>/`` 子目录 = 一个独立物性作业，交给通用引擎（cwd=子目录）。

        Note:
            potential/ 逐 epoch 不同，故 pair_style 的 hdnnp dir 用运行期 ``$(pwd)/potential``
            （引擎保证 cwd=本 epoch 目录），一条 cmd 即可服务所有 epoch。launcher 用 srun：
            本集群 LAMMPS 链接 Intel MPI（mpi/intel/17.0.7-thc），与 VASP 一致须用 srun；
            ``mpirun -np 1`` 时 Intel MPI 的 pmi_proxy 在共享节点上会间歇
            ``malloc(): memory corruption`` 崩溃（exit 255）。runner 读 ``$MY_LAUNCHER``
            （if_use_my_launcher=True 时引擎导出 ``srun -n $SLURM_NTASKS``）。

        Args:
            dir_run: n2p2 训练运行目录。
            if_sbatch: 是否真正提交 Slurm 作业。
            dict_args_to_submit: 传递给通用 Slurm 提交器的参数字典。

        Returns:
            物性测试根目录。

        Raises:
            FileNotFoundError: 缺少训练输入、scaling.data、权重文件或 LAMMPS 工具目录。
        """
        os.chdir(self.dir_root)
        dir_run = Path(dir_run)
        dir_post = dir_run.parent.parent / 'y_post' / dir_run.name
        dir_props = dir_post / 'properties'
        if dir_props.is_dir():
            shutil.rmtree(dir_props)

        os.makedirs(dir_props, exist_ok=True)

        # 1. 校验 + cutoff（pair_style hdnnp 近邻表 cutoff 必须 >= 最大对称函数 rc，留 0.01 Å 余量）
        src_input_nn = dir_run / 'input.nn'
        if not src_input_nn.is_file():
            raise FileNotFoundError(f"❌ Missing {src_input_nn}. Training not finished?")
        if not (dir_run / 'scaling.data').is_file():
            raise FileNotFoundError(f"❌ Missing {dir_run / 'scaling.data'}. nnp-train 应已生成 scaling.data。")
        if not self.dir_lmp_utils.is_dir():
            raise FileNotFoundError(f"❌ Missing lmp_utils source dir {self.dir_lmp_utils}.")
        cutoff = get_largest_rc_from_input_nn(str(src_input_nn)) + 0.01

        # 2. 收集 epoch（weights.<elem>.<epoch>.out -> <epoch> 去重排序）
        lepoch = sorted({w.name.split('.')[-2] for w in dir_run.glob('weights.*.*.out')})
        if not lepoch:
            raise FileNotFoundError(f"❌ No weights.*.*.out in {dir_run}. Training not finished?")

        # 3. 脚手架：建 y_epoch_scan/y_dir/<epoch>/（已存在则删除，遵守不递归删除约束）
        dir_epoch_scan_root = dir_props / 'y_epoch_scan'
        for epoch in lepoch:
            dir_epoch = dir_epoch_scan_root / 'y_dir' / epoch
            if dir_epoch.is_dir():
                shutil.rmtree(dir_epoch)
            lw = sorted(dir_run.glob(f'weights.*.{epoch}.out'))
            self._prepare_epoch_dir(dir_epoch, dir_run, lw, self.dir_lmp_utils)
        print(f"Prepared {len(lepoch)} epoch dir(s) under {dir_epoch_scan_root} "
              f"(hdnnp cutoff {cutoff:.2f}, pair_coeff {' '.join(self.lele)}).")


        # 4. 生成 pei_lmp_run_properties 的运行参数（脚本据此 sed 模板 + 写 general_mass.mod，cwd=本 epoch 目录）：
        #    $1 pair_style   整串 hdnnp：cutoff 已由 input.nn 算出；dir 用运行期 $(pwd)/potential（potential/ 逐 epoch 不同）
        #    $2 pair_coeff   元素串（hdnnp 用元素名，非势文件路径）
        #    $3 python_path  dft 环境（post/*.py 与 gsfe slab 建模需 mymetal）
        #    $4 mass_content general_mass.mod 内容：每元素一行 mass <type> <amu>，type 序与 pair_coeff 一致
        #    $5 lmp_template 本 epoch 的 ./template（= $(pwd)/template，用 runner 默认值，不在此传）
        from ase.data import atomic_numbers, atomic_masses
        python_path = self.path_python
        pair_style = (f'hdnnp {cutoff:.2f} dir $(pwd)/potential '
                      'showew no showewsum 100 resetew yes maxew 1000000 cflength 1.0 cfenergy 1.0')
        pair_coeff = ' '.join(self.lele)
        # mass_content：多元素用 \n 连接（runner 端 printf '%b' 展开成多行写入 general_mass.mod）
        mass_content = '\\n'.join(f'mass {i + 1} {atomic_masses[atomic_numbers[el]]:.4f}'
                                  for i, el in enumerate(self.lele))
        # on the PATH, not in cwd, so that it can be shared by all epochs
        runner = 'pei_lmp_run_properties'
        cmd = (f'{runner} "{pair_style}" "{pair_coeff}" '
               f'"{python_path}" "{mass_content}"')
        # 注入运行期参数（preset / 资源项由 dict_args_to_submit 提供，不在此固化）：
        #   cmd       runner + 4 个参数（见上；含 $(pwd)，必须留到作业里展开）
        #   path_root 扫描根（绝对）；其下 dir_root 的每个子目录 = 一个 epoch 物性作业
        #   dir_root  epoch 子目录的父目录（相对 path_root），覆盖 preset 的 ./y_dir
        # 用 {**d, ...} 生成新 dict（不改可变默认参数）；--if_sbatch 决定真提交还是 dry-run。
        args_to_submit = {**dict_args_to_submit,
                          'cmd': cmd,
                          'path_root': str(dir_epoch_scan_root.resolve()),
                          'dir_root': './y_dir',
                          'if_sbatch': True if if_sbatch else False,
                          }
        cli_args = self._dict_to_cli_args({**args_to_submit})


        self.last_jobids = []
        if if_sbatch:
            # pei_slurm_univ_submit 不回传作业号；each-subdir 的默认 shared 布局只提交
            # 1 个编排父作业，per-chunk / single-alloc 则可能提交多个。用提交前后的
            # squeue 差集捕获本次新增 Slurm 作业；只要 shared 父作业仍在等待，它就能
            # 覆盖后续尚未提交的 chunk 子作业尾部。
            before = self._snapshot_jobids()
            os.system(f'pei_slurm_univ_submit.py {cli_args}')
            after = self._snapshot_jobids()
            self.last_jobids = sorted(after - before)
            print(f"  post_properties observed {len(self.last_jobids)} new Slurm job(s): {self.last_jobids}")
        return dir_props


    def post_epoch_scan(self, dir_run: Path = Path('./train/y_n2p2_train/y_dir/001'),
                        dir_dft_root: Path = None) -> Path:
        """汇总跨 epoch 的 LAMMPS 物性扫描结果。

        把 ``properties/y_epoch_scan/y_dir/<epoch>/`` 各 epoch 的物性结果汇总到与 y_dir
        同级的 ``y_epoch_scan/`` 下，三类产物各自独立成表 + 图：

        - ``p_post_epoch_stretch.txt/pdf``：每 epoch 各相平衡 a/c、c/a(HCP)、E/atom、ΔE(相-FCC)。
        - ``p_post_epoch_cij.txt/pdf``：每 epoch 各相弹性常数 Cij。
        - ``p_post_epoch_gsfe.txt/pdf``：每 epoch 各滑移系的不稳定/稳定层错能（usf/sf）。

        物理量全部经 mymetal 的逐 epoch post 读取器取得（单一数据来源，避免重复解析逻辑）：
        ``my_read_stretch`` 读 p_post_stretch.txt 的 Extr infos(E)/rvector(a,c)，
        ``read_cij_energy`` 读 y_post_cij_energy.txt 的 C11..C44，
        ``read_output`` 从 y_post_gsfe.txt 的 gamma 表强制取 usf=max(gamma)、sf=最后一个 gamma。

        Note:
            某 epoch 缺 stretch/cij 的任一相则该表跳过此 epoch；gsfe 按滑移系逐个读取，
            缺失留 NaN。三表行数与 epoch 数不一致时只告警、不报错。

        Args:
            dir_run: n2p2 训练运行目录。
            dir_dft_root: DFT (VASP) 计算归档根目录（如 construct_dataset/calculate）。
                给定时用 :func:`read_dft_reference` 读出同样的物理量，叠成灰色虚线参考线
                画进三张 epoch 扫描图，并把读到的值打印到屏幕；为 None 时只画 LAMMPS 曲线。

        Returns:
            epoch 扫描汇总目录。

        Raises:
            FileNotFoundError: 缺少 properties/y_epoch_scan/y_dir 目录，或给定的 dir_dft_root 不存在。
        """

        from contextlib import redirect_stdout
        from io import StringIO
        from mymetal.post.stretch import my_read_stretch
        from mymetal.post.Cij_energy import read_cij_energy
        from mymetal.post.gsfe import read_output

        os.chdir(self.dir_root)
        dir_run = Path(dir_run)
        dir_props = dir_run.parent.parent / 'y_post' / dir_run.name / 'properties'
        dir_scan = dir_props / 'y_epoch_scan'      # 汇总表/图与 y_dir 同级（不放进 y_dir 内部）
        scan = dir_scan / 'y_dir'
        if not scan.is_dir():
            raise FileNotFoundError(f"❌ Missing {scan}. Run post_properties(if_sbatch=True) first.")

        lats, cij_keys = _LATS, _CIJ_KEYS
        gsfe_types, all_types = _GSFE_TYPES, _ALL_GSFE_TYPES

        epochs = sorted([int(d.name) for d in scan.iterdir() if d.is_dir() and d.name.isdigit()])

        rows_stretch, rows_cij, rows_gsfe = [], [], []
        for ep in epochs:
            epp = f'{ep:06d}'

            # 1) stretch：各相平衡 a/c (Extr rvector) 与 E/atom (Extr infos)，缺任一相则跳过本 epoch
            rec, ok = {'epoch': ep}, True
            for lat in lats:
                ps = scan / epp / 'y_stretch' / lat / 'p_post_stretch.txt'
                if not ps.is_file():
                    ok = False; break
                with redirect_stdout(StringIO()):              # my_read_stretch 会回显注释头，抑制
                    *_, (_, _, extr_y, extr_rvec) = my_read_stretch(str(ps))
                rec[f'a_{lat}'] = float(extr_rvec[0])
                rec[f'c_{lat}'] = float(extr_rvec[2])
                rec[f'E_{lat}'] = float(extr_y)
            if ok:
                rec['ca_hcp'] = rec['c_hcp'] / rec['a_hcp']
                rec['dE_hcp_fcc'] = (rec['E_hcp'] - rec['E_fcc']) * 1e3
                rec['dE_bcc_fcc'] = (rec['E_bcc'] - rec['E_fcc']) * 1e3
                rows_stretch.append(rec)
            else:
                print(f'skip epoch {epp}: incomplete stretch')

            # 2) cij：各相 C11..C44，缺任一相则跳过本 epoch
            rec, ok = {'epoch': ep}, True
            for lat in lats:
                pc = scan / epp / 'y_Cij_energy' / lat / 'y_post_cij_energy.txt'
                if not pc.is_file():
                    ok = False; break
                cij = read_cij_energy(str(pc))
                for k in cij_keys:
                    rec[f'{k}_{lat}'] = cij[k]
            if ok:
                rows_cij.append(rec)
            else:
                print(f'skip epoch {epp}: incomplete cij')

            # 3) gsfe：逐滑移系读取 usf=max(gamma)、sf=最后一个 gamma，缺失留 NaN
            rec, ok = {'epoch': ep}, True
            for phase, ltype in gsfe_types.items():
                for t in ltype:
                    pg = scan / epp / 'y_gsfe' / phase / t / 'y_post_gsfe.txt'
                    if not pg.is_file():
                        ok = False; break
                    g = read_output(str(pg))
                    rec[f'usf_{t}'] = g['usf_max'] if g['usf_max'] is not None else np.nan
                    rec[f'sf_{t}'] = g['sf_min'] if g['sf_min'] is not None else np.nan
            if ok:
                rows_gsfe.append(rec)
            else:
                print(f'skip epoch {epp}: incomplete gsfe')

        if len(rows_stretch) != len(epochs) or len(rows_cij) != len(epochs) or len(rows_gsfe) != len(epochs):
            print(f"⚠️ Warning: Inconsistent epoch counts: "
                             f"stretch {len(rows_stretch)}, cij {len(rows_cij)}, gsfe {len(rows_gsfe)}, epochs {len(epochs)}. "
                             f"Check {scan} for missing files.")

        # 列顺序（缺的列在 reindex 时自动补 NaN，不影响表/图）
        cols_stretch, cols_cij, cols_gsfe = _COLS_STRETCH, _COLS_CIJ, _COLS_GSFE

        nstretch, n_cij, n_gsfe = 0, 0, 0

        # DFT 参考基线（construct_dataset 归档）：与逐 epoch 量同读取器、同键名，叠成灰色虚线参考。
        # 给定 dir_dft_root 才读；目录不存在直接报错（结构性前提），缺单个文件由 read_dft_reference 跳过。
        dft = None
        if dir_dft_root is not None:
            if not Path(dir_dft_root).is_dir():
                raise FileNotFoundError(f"❌ Missing DFT archive root {dir_dft_root}.")
            dft = read_dft_reference(dir_dft_root)
            print(f"================ 📊 DFT reference ({dir_dft_root})")
            for group in ('stretch', 'cij', 'gsfe'):
                print(f"  [{group}] {len(dft[group])} value(s)")
                for k, v in dft[group].items():
                    print(f"    {k:16s} = {v:.6g}")

        if rows_stretch:
            df_stretch = pd.DataFrame(rows_stretch).sort_values('epoch').reset_index(drop=True)[cols_stretch]
            self._write_table(dir_scan / 'p_post_epoch_stretch.txt',
                            ['# LAMMPS (pair_style hdnnp) equilibrium stretch properties vs training epoch',
                            '# FCC/BCC a/c are conventional cubic lengths; HCP a/c are hexagonal; E (eV/atom), dE (meV/atom)'],
                            df_stretch, float_format='14.6f')
            my_plot_epoch_stretch(df_stretch, dir_scan / 'p_post_epoch_stretch.pdf',
                                  dft=dft['stretch'] if dft else None)
            nstretch = len(df_stretch)
        else:
            print('No complete stretch epochs; skip stretch.txt/pdf.')

        if rows_cij:
            df_cij = pd.DataFrame(rows_cij).sort_values('epoch').reset_index(drop=True)[cols_cij]
            self._write_table(dir_scan / 'p_post_epoch_cij.txt',
                              ['# LAMMPS (pair_style hdnnp) elastic constants Cij vs training epoch',
                               '# Cij (GPa); cubic phases use C11=C33, C12=C13'],
                              df_cij, float_format='12.4f')
            my_plot_epoch_cij(df_cij, dir_scan / 'p_post_epoch_cij.pdf',
                              dft=dft['cij'] if dft else None)
            n_cij = len(df_cij)
        else:
            print('No complete cij epochs; skip cij.txt/pdf.')

        if rows_gsfe:
            df_gsfe = pd.DataFrame(rows_gsfe).sort_values('epoch').reset_index(drop=True).reindex(columns=cols_gsfe)
            self._write_table(dir_scan / 'p_post_epoch_gsfe.txt',
                              ['# LAMMPS (pair_style hdnnp) stacking-fault energies vs training epoch',
                               '# usf = max(gamma), sf = last gamma value; units mJ/m^2'],
                              df_gsfe, float_format='12.4f')
            my_plot_epoch_gsfe(df_gsfe, dir_scan / 'p_post_epoch_gsfe.pdf', types=all_types,
                               dft=dft['gsfe'] if dft else None)
            n_gsfe = len(df_gsfe)
        else:
            print('No gsfe epochs; skip gsfe.txt/pdf.')

        print(f"✅ Aggregated epochs -> {dir_scan}: "
              f"stretch={nstretch}, cij={n_cij}, gsfe={n_gsfe} "
              f"(stretch/cij/gsfe .txt + .pdf, sibling of y_dir)")
        return dir_scan


    def check_interface(self, dir_run: Path = Path('./train/y_n2p2_train/y_dir/001'),
                        epoch: int = None, if_run: bool = False) -> Path:
        """检查 nnp-predict 与 LAMMPS hdnnp 的单点预测一致性（接口一致性检查 C0）。

        对同一结构分别用 nnp-predict（训练工具链，n2p2 2.3.0）和 LAMMPS
        ``pair_style hdnnp``（LAMMPS 内嵌 2.2.0）算单点，比较总能量（目标 ~1e-5 eV）和
        逐原子力分量（目标 ~1e-5 eV/Å），确认两套 n2p2 对同一 nnp-data 给出一致预测——
        力还额外覆盖 LAMMPS 的力累加 + ghost 原子通信路径（能量检查压不到）。产物落在
        ``y_post/<run>/check_interface/``：nnp/（nnp-predict 工作目录）、lmp/（LAMMPS）、
        run_check.bash、p_post_check_interface.txt。

        Note:
            刻意取一个轴对齐正交胞：能量与旋转无关，且 ASE 写 lammps-data 不旋转坐标系
            -> 两边力同系、逐分量直接比，无需坐标变换。``if_run=True`` 时在登录节点直接跑
            两边（都是单点、秒级）并写结论；否则只准备，由用户手动跑 run_check.bash。

        Args:
            dir_run: n2p2 训练运行目录。
            epoch: 用于检查的 epoch；为 None 时读取 post_training 写出的 epoch。
            if_run: 是否直接运行接口检查脚本并解析结果。

        Returns:
            接口检查目录。

        Raises:
            FileNotFoundError: 缺少 epoch 信息、权重文件、势函数文件或检查脚本。
        """
        import numpy as _np
        from ase.io import write as _ase_write

        os.chdir(self.dir_root)
        dir_run = Path(dir_run)
        dir_post = dir_run.parent.parent / 'y_post' / dir_run.name
        dir_chk = dir_post / 'check_interface'
        dnnp, dlmp = dir_chk / 'nnp', dir_chk / 'lmp'
        os.makedirs(dnnp, exist_ok=True)
        os.makedirs(dlmp, exist_ok=True)

        if epoch is None:
            pe = dir_post / 'training' / 'epoch.txt'
            if not pe.is_file():
                raise FileNotFoundError(f"❌ Missing {pe}. Run post_training first or pass epoch=.")
            epoch = int(pe.read_text(encoding='utf-8').strip())

        # 部署势函数三件套到 nnp/（nnp-predict 按 cwd 读；LAMMPS hdnnp 的 dir 也指向它）
        lw = sorted(dir_run.glob(f'weights.*.{epoch:06d}.out'))
        if not lw:
            raise FileNotFoundError(f"❌ No weights.*.{epoch:06d}.out in {dir_run}.")
        for fn in ['input.nn', 'scaling.data']:
            shutil.copy(dir_run / fn, dnnp / fn)
        for w in lw:
            shutil.copy(w, dnnp / f"weights.{w.name.split('.')[1]}.data")

        # 取第一个正交胞、原子数>=4 的结构（fcc/bcc 常规胞）
        # 正交不是物理必需，而是消除三斜→LAMMPS 几何转换这个混淆因素，让接口一致性检查最干净。
        d = nnpdata()
        d.load_from_datafile(str(dir_run / 'input.data'))

        def _is_orth(c, atol=1e-6):
            c = _np.array(c)
            return _np.max(_np.abs(c - _np.diag(_np.diag(c)))) < atol

        idx = 0
        for i, a in enumerate(d.latoms):
            if _is_orth(a.cell) and len(a) >= 4:
                idx = i
                break
        atoms = d.latoms[idx]
        e_dft = atoms.get_potential_energy()
        sym = atoms.get_chemical_symbols()
        elements = list(dict.fromkeys(sym))   # 去重并保留首次出现顺序 -> LAMMPS type 顺序

        # nnp/input.data（单结构，Å/eV）：直接用 mymetal nnpdata 的写出器，
        # 与训练集 input.data 同一套格式/精度（lattice、atom+force、energy、charge）。
        d.write_from_ase(str(dnnp / 'input.data'), atoms, tag=d.ltags[idx],
                         file_name=f'check_interface_src_index_{idx}', append=False)

        # lmp/lmp.data（ASE 把三斜胞转 LAMMPS 约定；能量不变）+ lmp.in
        # masses=True 让 ASE 写出 Masses 段（按元素质量），lmp.in 不再硬编码 mass；
        # specorder / pair_coeff 由结构元素自适应得到，与 data 文件的 type 顺序一致。
        _ase_write(str(dlmp / 'lmp.data'), atoms, format='lammps-data',
                   specorder=elements, masses=True)
        dnnp_abs = str(dnnp.resolve())
        # pair_style hdnnp 近邻表 cutoff 必须 >= 最大对称函数 rc；从 input.nn 读最大 rc 再留 0.01 Å 余量。
        rc_max = get_largest_rc_from_input_nn(str(dnnp / 'input.nn'))
        cutoff = rc_max + 0.01
        in_c0 = "\n".join([
            "units           metal",
            "atom_style      atomic",
            "boundary        p p p",
            "read_data       lmp.data",
            f"pair_style      hdnnp {cutoff:.2f} dir {dnnp_abs} showew yes showewsum 1 resetew yes maxew 1000000 cflength 1.0 cfenergy 1.0",
            f"pair_coeff      * * {' '.join(elements)}",
            "dump            f all custom 1 forces.dump id fx fy fz",
            'dump_modify     f sort id format line "%d %.16e %.16e %.16e"',
            "run 0",
            "variable e equal pe",
            'print "LAMMPS_PE_eV $(v_e:%.12f)"',
            "",
        ])
        (dlmp / 'lmp.in').write_text(in_c0, encoding='utf-8')

        # run_check.bash 与本机环境强耦合（模块名、二进制路径），作为脚本放在 file/，这里只部署。
        src_run_check = self.dir_file / 'run_check.bash'
        if not src_run_check.is_file():
            raise FileNotFoundError(f"❌ Missing {src_run_check} (env-specific interface-check runner).")
        shutil.copy(src_run_check, dir_chk / 'run_check.bash')

        print(f"Prepared interface check in {dir_chk}: structure index {idx} "
              f"(tag {d.ltags[idx]}, {len(atoms)} atoms), epoch {epoch:06d}.")

        if if_run:
            os.system(f'bash {dir_chk / "run_check.bash"}')
            e_nnp = e_lmp = None
            od = dnnp / 'output.data'
            if od.is_file():
                for ln in od.read_text(encoding='utf-8').splitlines():
                    if ln.strip().startswith('energy'):
                        e_nnp = float(ln.split()[1])
            import re as _re
            ls = dlmp / 'lmp.stdout'
            if ls.is_file():
                mt = _re.search(r'LAMMPS_PE_eV\s+(-?\d+\.\d+)', ls.read_text(encoding='utf-8'))
                if mt:
                    e_lmp = float(mt.group(1))

            # 力：nnp-predict 的 nnforces.out 前 3 列（NNP 力，按原子序）vs LAMMPS forces.dump（sort id）。
            # 已选轴对齐正交胞 -> ASE 不旋转坐标系，两边同系，逐分量直接比，无需任何坐标变换。
            def _read_nnforces(p):
                if not p.is_file():
                    return None
                rows = [ln.split()[:3] for ln in p.read_text(encoding='utf-8').splitlines()
                        if ln.strip() and not ln.lstrip().startswith('#')]
                return _np.array(rows, dtype=float) if rows else None

            def _read_lmp_dump(p):
                if not p.is_file():
                    return None
                lines = p.read_text(encoding='utf-8').splitlines()
                for i, ln in enumerate(lines):
                    if ln.startswith('ITEM: ATOMS'):     # 列：id fx fy fz（dump_modify sort id）
                        rows = []
                        for l in lines[i + 1:]:
                            if l.startswith('ITEM:'):     # 只取第一帧
                                break
                            if l.strip():
                                rows.append(l.split())
                        return _np.array([[float(c) for c in r[1:4]] for r in rows]) if rows else None
                return None

            F_nnp = _read_nnforces(dnnp / 'nnforces.out')
            F_lmp = _read_lmp_dump(dlmp / 'forces.dump')

            with open(dir_chk / 'p_post_check_interface.txt', 'w', encoding='utf-8') as f:
                f.write('# Interface consistency check: nnp-predict (2.3.0) vs LAMMPS hdnnp (2.2.0)\n')
                f.write(f'structure   src_index={idx}, tag={d.ltags[idx]}, natoms={len(atoms)}, epoch={epoch:06d}\n')
                if e_nnp is not None and e_lmp is not None:
                    dlt = abs(e_nnp - e_lmp)
                    f.write(f'E_nnp_predict  {e_nnp:.12f} eV\n')
                    f.write(f'E_LAMMPS_hdnnp {e_lmp:.12f} eV\n')
                    f.write(f'|delta_E|      {dlt:.3e} eV  ({dlt/len(atoms)*1e3:.3e} meV/atom)\n')
                    f.write(f'E_DFT_ref      {e_dft:.12f} eV\n')
                    f.write(f'verdict_E      {"✅ PASS" if dlt < 1e-4 else "❌ FAIL"} (target ~1e-5 eV)\n')
                    print(f"  nnp-predict {e_nnp:.9f} eV  vs  LAMMPS {e_lmp:.9f} eV  |delta|={dlt:.2e} eV "
                          f"-> {'✅ PASS' if dlt < 1e-4 else '❌ FAIL'}")
                else:
                    f.write('ERROR: could not parse one/both energies; check nnp/ and lmp/ logs.\n')
                    print("  ❌ could not parse energies; check nnp/nnp-predict.stdout and lmp/lmp.stdout")

                # 力比较（逐分量；nnforces.out 为 %16.8E -> 比较精度地板 ~1e-7，阈值取 ~1e-5 eV/A）
                if F_nnp is not None and F_lmp is not None and F_nnp.shape == F_lmp.shape:
                    maxdF = float(_np.abs(F_nnp - F_lmp).max())
                    maxF = float(_np.abs(F_nnp).max())
                    f.write(f'max|F|         {maxF:.6e} eV/A  (nnp-predict, force scale)\n')
                    f.write(f'max|dF_comp|   {maxdF:.3e} eV/A  (component-wise, orthogonal cell, no transform)\n')
                    f.write(f'verdict_F      {"✅ PASS" if maxdF < 1e-5 else "❌ FAIL"} (target <~1e-5 eV/A)\n')
                    print(f"  forces: max|dF|={maxdF:.2e} eV/A (max|F|={maxF:.2e}) "
                          f"-> {'✅ PASS' if maxdF < 1e-5 else '❌ FAIL'}")
                else:
                    sh_n = None if F_nnp is None else F_nnp.shape
                    sh_l = None if F_lmp is None else F_lmp.shape
                    f.write(f'ERROR: could not parse/match forces (nnp {sh_n}, lmp {sh_l}); '
                            'check nnp/nnforces.out and lmp/forces.dump.\n')
                    print(f"  ❌ could not compare forces (nnp {sh_n}, lmp {sh_l})")
        return dir_chk


    # summary y_post 下的所有信息：check_interface、properties、training
    @staticmethod
    def _read_epoch_table(path: Path) -> pd.DataFrame:
        """读回 ``post_epoch_scan`` 写出的 ``p_post_epoch_*.txt`` 定宽表。

        ``_write_table`` 的格式是若干 ``#`` 注释行 + 列名行 + 右对齐的定宽数据行，
        故按空白切分、丢弃注释即可还原 DataFrame（含 ``epoch`` 列）。

        Args:
            path: p_post_epoch_{stretch,cij,gsfe}.txt 之一。

        Returns:
            以原列名为列的 DataFrame。
        """
        return pd.read_csv(path, sep=r'\s+', comment='#')


    @staticmethod
    def _ensemble_stats(ldf: list, cols: list, index_col: str = 'epoch',
                        band: str = 'minmax') -> tuple:
        """把多个 run 的同一张表叠成 across-run 的均值 / 包络 / 标准差。

        各 run 的索引键（epoch 或 tag）取并集后对齐（某 run 缺的行或列补 NaN），再沿
        run 轴做 nan-aware 统计，因此 run 之间行数不一致（例如某 run 少了 epoch 0）
        不会丢掉其余 run 在该行上的信息。

        Args:
            ldf: 各 run 的 DataFrame 列表，每个含 ``index_col`` 列。
            cols: 期望的列顺序（可含 ``index_col``，会被剔除）；缺的列补 NaN。
            index_col: 对齐用的索引列名，``'epoch'``（epoch 扫描 / 学习曲线）或
                ``'tag'``（按 tag 的 RMSE）。
            band: 色带取法。``'minmax'`` 用逐行的跨 run 上下限（默认，对应 chap_4
                图 1.3/1.5 中的“上下限”读法）；``'std'`` 用 mean ± std。

        Returns:
            tuple: ``(df_mean, df_lo, df_hi, df_out)``。前三者以 ``index_col`` 打头、
            列名与 ``cols`` 一致，直接喂给绘图函数；``df_out`` 是写盘用的宽表
            （``index_col``、``n_runs``，以及每个量的 ``_mean/_min/_max/_std``）。

        Raises:
            ValueError: band 不是 'minmax' 或 'std'。
        """
        if band not in _BAND_LABEL:
            raise ValueError(f"band must be one of {sorted(_BAND_LABEL)}, got {band!r}")

        lcol = [c for c in cols if c != index_col]
        lkey = sorted(set().union(*(set(df[index_col]) for df in ldf)))

        # (n_run, n_key, n_col)：对齐到索引并集 x 期望列，缺失处为 NaN
        arr = np.full((len(ldf), len(lkey), len(lcol)), np.nan)
        for i, df in enumerate(ldf):
            arr[i] = df.set_index(index_col).reindex(index=lkey, columns=lcol).to_numpy(dtype=float)

        # 全 NaN 的 (key, col) 切片会让 nanmean/nanmin 报 RuntimeWarning 并返回 NaN；
        # NaN 正是我们想要的“该点无数据”，故只静音警告。std 用 ddof=1（样本标准差），
        # 单 run 时退化为 NaN —— 与“无法估计离散度”一致。
        with np.errstate(invalid='ignore'), warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            mean = np.nanmean(arr, axis=0)
            vmin = np.nanmin(arr, axis=0)
            vmax = np.nanmax(arr, axis=0)
            std = np.nanstd(arr, axis=0, ddof=1) if len(ldf) > 1 else np.full_like(mean, np.nan)

        lo, hi = (vmin, vmax) if band == 'minmax' else (mean - std, mean + std)
        # n_runs：该行有几个 run 提供了数据（整行至少一个有限值）
        n_runs = np.isfinite(arr).any(axis=2).sum(axis=0)

        def _frame(a):
            return pd.DataFrame(a, columns=lcol).assign(**{index_col: lkey})[[index_col] + lcol]

        df_out = pd.DataFrame({index_col: lkey, 'n_runs': n_runs})
        for j, c in enumerate(lcol):                        # 每个量 4 列，保持 cols 的顺序
            df_out[f'{c}_mean'] = mean[:, j]
            df_out[f'{c}_min'] = vmin[:, j]
            df_out[f'{c}_max'] = vmax[:, j]
            df_out[f'{c}_std'] = std[:, j]

        return _frame(mean), _frame(lo), _frame(hi), df_out


    @staticmethod
    def _read_check_interface(path: Path) -> dict:
        """从 ``p_post_check_interface.txt`` 抽出 verdict 与量化偏差。

        ``check_interface`` 写的是人读格式，不是表格；本函数按行首关键字取值：
        ``verdict_E`` / ``verdict_F`` 行含 ``PASS`` 或 ``FAIL``（解析失败时该行整个
        缺失，此时留 NaN，绘图上表现为空缺而不是默认通过）；``|delta_E|`` 与
        ``max|dF_comp|`` 行的第二个字段是数值。

        Args:
            path: 某个 run 的 p_post_check_interface.txt。

        Returns:
            dict: ``epoch``、``tag``、``verdict_E``/``verdict_F``（1=PASS，2=FAIL，
            缺失为 NaN）、``delta_E_eV``、``max_dF_comp_eV_A``。
        """
        rec = {'epoch': np.nan, 'tag': '', 'verdict_E': np.nan, 'verdict_F': np.nan,
               'delta_E_eV': np.nan, 'max_dF_comp_eV_A': np.nan}
        for line in Path(path).read_text(encoding='utf-8').splitlines():
            lfield = line.split()
            if not lfield:
                continue
            head = lfield[0]
            if head == 'structure':                        # structure src_index=.., tag=.., epoch=..
                for kv in line.replace(',', ' ').split():
                    if kv.startswith('tag='):
                        rec['tag'] = kv.split('=', 1)[1]
                    elif kv.startswith('epoch='):
                        rec['epoch'] = int(kv.split('=', 1)[1])
            elif head in ('verdict_E', 'verdict_F'):
                if 'PASS' in line:
                    rec[head] = VERDICT_CODE['PASS']
                elif 'FAIL' in line:
                    rec[head] = VERDICT_CODE['FAIL']
            elif head == '|delta_E|':
                rec['delta_E_eV'] = float(lfield[1])
            elif head == 'max|dF_comp|':
                rec['max_dF_comp_eV_A'] = float(lfield[1])
        return rec


    def _summary_root(self, ldir_run: list, band: str = None) -> tuple:
        """校验汇总入参并解析出 run 目录列表与写盘根目录（各 run 共同的 ``y_post/``）。

        三个 ``post_*_summary`` 共用：先校验再落盘，参数写错时不产生任何文件。

        Args:
            ldir_run: 参与汇总的 n2p2 训练运行目录列表。
            band: 若给定则一并校验其取值合法。

        Returns:
            tuple: ``(ldir_run, dir_post_root)``，均为 Path。

        Raises:
            ValueError: 未给 ldir_run，或 band 不是 'minmax'/'std'。
        """
        os.chdir(self.dir_root)
        if not ldir_run:
            raise ValueError("ldir_run must be a non-empty list of training run dirs.")
        if band is not None and band not in _BAND_LABEL:
            raise ValueError(f"band must be one of {sorted(_BAND_LABEL)}, got {band!r}")
        ldir_run = [Path(d) for d in ldir_run]
        return ldir_run, ldir_run[0].parent.parent / 'y_post'


    def _collect_run_tables(self, ldir_run: list, dir_post_root: Path, rel_path: Path) -> tuple:
        """读取每个 run 的 ``y_post/<run-id>/<rel_path>`` 表，缺文件的 run 跳过并告警。

        Args:
            ldir_run: 训练运行目录列表（只用其 ``name`` 定位 y_post 子目录）。
            dir_post_root: ``y_post/`` 根目录。
            rel_path: 相对 ``y_post/<run-id>/`` 的表文件路径。

        Returns:
            tuple: ``(ldf, lrid)`` —— 成功读到的 DataFrame 列表与对应的 run-id 列表。
        """
        ldf, lrid = [], []
        for dir_run in ldir_run:
            p = dir_post_root / dir_run.name / rel_path
            if not p.is_file():
                print(f'skip run {dir_run.name}: missing {p}')
                continue
            ldf.append(self._read_epoch_table(p))
            lrid.append(dir_run.name)
        return ldf, lrid


    @staticmethod
    def _scan_epoch_absence(ldir_run: list, dtable: dict) -> pd.DataFrame:
        """逐 (run, kind, epoch) 找出 epoch 扫描三张表里缺了什么。

        epoch 全集取**跨 run 跨 kind 的并集**：只有某个 run 的某张表独有的 epoch，才能
        暴露出其它 run / 其它 kind 在该 epoch 上的缺失。这与 :meth:`_ensemble_stats` 的
        对齐口径一致，故本表列出的每一条都对应汇总里一个 NaN 参与统计的点。

        Args:
            ldir_run: 参与汇总的 run 目录列表（只用其 ``name``）。
            dtable: ``kind -> {run-id: DataFrame}``，post_epoch_scan_summary 已读到的表。

        Returns:
            以 :data:`_COLS_ABSENCE` 为列的 DataFrame，无缺失时为空表。
        """
        leps = sorted(set().union(*(set(df['epoch']) for ddf in dtable.values() for df in ddf.values())))

        rows = []
        for dir_run in ldir_run:
            rid = dir_run.name
            for kind, (cols, *_) in _EPOCH_KINDS.items():
                df = dtable.get(kind, {}).get(rid)
                if df is None:
                    rows.append({'run': rid, 'kind': kind, 'epoch': -1,
                                 'status': 'missing_table', 'item': f'p_post_epoch_{kind}.txt'})
                    continue

                lcol = [c for c in cols if c != 'epoch']
                # reindex 到 epoch 全集 x 期望列：缺行、缺列、写进去的 NaN 统一变成 NaN
                d = df.set_index('epoch').reindex(index=leps, columns=lcol).apply(pd.to_numeric, errors='coerce')
                bad = ~np.isfinite(d.to_numpy(dtype=float))
                have = set(df['epoch'])
                for i, ep in enumerate(leps):
                    if ep not in have:                     # 整行没进表，逐列展开只会刷屏
                        rows.append({'run': rid, 'kind': kind, 'epoch': ep,
                                     'status': 'missing_epoch', 'item': '*'})
                        continue
                    for j in np.flatnonzero(bad[i]):
                        rows.append({'run': rid, 'kind': kind, 'epoch': ep,
                                     'status': 'missing_value', 'item': lcol[j]})
        return pd.DataFrame(rows, columns=_COLS_ABSENCE)


    def _write_absence(self, path: Path, df: pd.DataFrame, ldir_run: list) -> None:
        """把 :meth:`_scan_epoch_absence` 的结果写成 ``p_post_absence.txt``。

        无缺失时也写文件（只有注释头），这样"文件不在"永远意味着汇总没跑过，而不是
        "跑过且什么都不缺"。

        Args:
            path: 输出路径（``y_post/p_post_absence.txt``）。
            df: 缺失记录表，可以为空。
            ldir_run: 参与汇总的 run 目录列表，写进表头备查。
        """
        header = ['# Missing / incomplete inputs behind the ensemble epoch-scan summary',
                  '# One row per (run, kind, epoch, item) that post_epoch_scan_summary had to treat as NaN',
                  f'# {len(ldir_run)} run(s): {", ".join(d.name for d in ldir_run)}',
                  '# status:'] + [f'#   {k:14s} {v}' for k, v in _ABSENCE_STATUS.items()]
        if df.empty:
            Path(path).write_text('\n'.join(header + ['# (none: every run has every epoch of every kind)', '']),
                                  encoding='utf-8')
        else:
            self._write_table(path, header, df)


    def post_epoch_scan_summary(self, ldir_run: list = None, dir_dft_root: Path = None,
                                dir_summary: Path = None, band: str = 'minmax') -> Path:
        """把多个 run 的 epoch 扫描曲线汇总成 ensemble 均值 + 上下限图（chap_4 图 1.3/1.5 风格）。

        :meth:`post_epoch_scan` 已为每个 run 在 ``y_post/<run-id>/properties/y_epoch_scan/``
        下写好 ``p_post_epoch_{stretch,cij,gsfe}.txt``；本方法只**读**这些表，跨 run 对齐
        epoch 后统计，把汇总产物直接写在 ``y_post/`` 下（与各 run 目录平级）：

        - ``p_post_epoch_stretch_summary.txt/pdf``
        - ``p_post_epoch_cij_summary.txt/pdf``
        - ``p_post_epoch_gsfe_summary.txt/pdf``
        - ``p_post_absence.txt``：缺失清单，逐 (run, kind, epoch, item) 记录本次汇总里
          被当成 NaN 的点 —— 缺整张表、缺某个 epoch 的行、或某一列是 NaN。无缺失时也写，
          只有注释头。

        图沿用逐 run 的面板布局（stretch 3x3、cij 5x3、gsfe 2x6），同一物理量在汇总图
        中的行列位置与它在单 run 图中的位置一致；每个面板画跨 run 均值（C0 线 + C1 圆点）
        并把上下限填成同色半透明色带，DFT 参考仍是灰色水平虚线。

        Note:
            纯读取 + 汇总，不触碰任何 run 自己的 ``y_post/<run-id>/`` 产物，可安全重跑。
            各 run 的 epoch 数允许不同（并集对齐、nan-aware 统计）。DFT 的 Cij 参考来自
            ``read_dft_reference`` 的默认 ``cij_subdir='y_cij_energy_small'``（小应变谐性
            弹性常数），与逐 epoch 的 LAMMPS 评估口径一致。

        Args:
            ldir_run: 参与汇总的 n2p2 训练运行目录列表（如 ``[.../y_dir/001, ...]``）；
                汇总输出落在这些 run 共同的 ``y_post/`` 下。
            dir_dft_root: DFT (VASP) 计算归档根目录；给定时叠加灰色虚线 DFT 参考。
            dir_summary: 汇总输出目录；默认就是 ``y_post/`` 本身。
            band: 色带取法，``'minmax'``（默认，逐 epoch 跨 run 上下限）或 ``'std'``
                （mean ± 样本标准差，对应 chap_4 图 1.3 的透明背景读法）。

        Returns:
            汇总输出目录。缺表的 run 只跳过并告警；三类全缺时只写 ``p_post_absence.txt``。

        Raises:
            ValueError: 未给 ldir_run，或 band 不是 'minmax'/'std'。
            FileNotFoundError: 给定的 dir_dft_root 不存在。
        """
        # 先校验入参再读 DFT / 落盘：参数写错时不产生任何文件
        ldir_run, dir_post_root = self._summary_root(ldir_run, band)
        dir_summary = Path(dir_summary) if dir_summary else dir_post_root

        # DFT 参考基线：与逐 run 图同一读取器、同一键名，故可直接叠成灰色水平虚线
        dft = None
        if dir_dft_root is not None:
            if not Path(dir_dft_root).is_dir():
                raise FileNotFoundError(f"❌ Missing DFT archive root {dir_dft_root}.")
            dft = read_dft_reference(dir_dft_root)
            print(f"================ 📊 DFT reference ({dir_dft_root})")
            for group in ('stretch', 'cij', 'gsfe'):
                print(f"  [{group}] {len(dft[group])} value(s)")
                for k, v in dft[group].items():
                    print(f"    {k:16s} = {v:.6g}")

        os.makedirs(dir_summary, exist_ok=True)            # 只保证存在；不动已算好的 y_post/<run-id>/

        lplot = {'stretch': my_plot_epoch_stretch, 'cij': my_plot_epoch_cij, 'gsfe': my_plot_epoch_gsfe}
        band_label = _BAND_LABEL[band]

        nwritten, dtable = 0, {}
        for kind, (cols, float_format, unit_note) in _EPOCH_KINDS.items():
            ldf, lrid = self._collect_run_tables(
                ldir_run, dir_post_root, Path('properties') / 'y_epoch_scan' / f'p_post_epoch_{kind}.txt')
            dtable[kind] = dict(zip(lrid, ldf))         # 留给 _scan_epoch_absence，缺表的 run 自然不在其中
            if not ldf:
                print(f'No run has p_post_epoch_{kind}.txt; skip {kind} summary.')
                continue

            lnep = {rid: len(df) for rid, df in zip(lrid, ldf)}
            if len(set(lnep.values())) > 1:                 # epoch 数不齐只是提示；统计已按并集对齐
                print(f"⚠️ {kind}: runs disagree on epoch count {lnep}; "
                      f"aligning on the union and using nan-aware statistics.")

            df_mean, df_lo, df_hi, df_out = self._ensemble_stats(ldf, cols, 'epoch', band)
            self._write_table(dir_summary / f'p_post_epoch_{kind}_summary.txt',
                              [f'# Ensemble summary of LAMMPS (pair_style hdnnp) {kind} vs training epoch',
                               f'# {unit_note}',
                               f'# {len(ldf)} run(s): {", ".join(lrid)}',
                               '# n_runs = runs contributing to this epoch; per quantity: _mean/_min/_max/_std '
                               '(std is the sample std, ddof=1, over runs)'],
                              df_out, float_format=float_format)

            kwargs = {'types': _ALL_GSFE_TYPES} if kind == 'gsfe' else {}
            lplot[kind](df_mean, dir_summary / f'p_post_epoch_{kind}_summary.pdf',
                        dft=dft[kind] if dft else None,
                        df_lo=df_lo, df_hi=df_hi, band_label=band_label.format(n=len(ldf)),
                        **kwargs)
            nwritten += 1
            print(f'  {kind}: {len(ldf)} run(s), {len(df_out)} epoch(s) -> '
                  f'p_post_epoch_{kind}_summary.txt/pdf')

        # 缺失清单：三类产物看完再写，才能用跨 kind 的 epoch 并集判定缺行
        df_absence = self._scan_epoch_absence(ldir_run, dtable)
        self._write_absence(dir_summary / 'p_post_absence.txt', df_absence, ldir_run)
        print(f'  absence: {len(df_absence)} record(s) -> p_post_absence.txt')

        print(f"✅ Ensemble epoch-scan summary ({band}) -> {dir_summary}: {nwritten}/3 group(s) written")
        return dir_summary

    def post_training_summary(self, ldir_run: list = None, dir_summary: Path = None,
                              band: str = 'minmax') -> Path:
        """跨 run 汇总 ``y_post/<run-id>/training/`` 的训练误差，写在 ``y_post/`` 下。

        两组产物，都用与 epoch 物性汇总相同的 ensemble 读法（均值 + 上下限）：

        - ``p_post_rmse_summary.txt/pdf``：训练集 E/F RMSE 随 epoch 的曲线，2 个 panel
          （上能量、下力）。数据来自各 run 的 ``p_post_learning_curve.txt`` —— 那是
          唯一含逐 epoch RMSE 的表；``p_post_rmse.txt`` 只有末 epoch 的 4 个标量，
          撑不起“均值线 + 色带 vs epoch”的画法（其数值等于本汇总里 by-tag 的 TOTAL 行）。
        - ``p_post_rmse_by_tag_summary.txt/pdf``：按 tag 的 E/F RMSE 柱状图，柱高为跨 run
          均值，误差棒为跨 run 上下限（或 mean ± std）。

        Note:
            纯读取 + 汇总，不触碰任何 run 自己的 ``y_post/<run-id>/`` 产物，可安全重跑。
            by-tag 表按 ``tag`` 对齐（不是按行号），故某 run 少一个 tag 也不会错位；
            ``n_struct`` / ``n_fcomp`` 是数据集属性而非训练结果，取第一个 run 的值透传，
            各 run 不一致时告警。

        Args:
            ldir_run: 参与汇总的 n2p2 训练运行目录列表。
            dir_summary: 汇总输出目录；默认就是 ``y_post/`` 本身。
            band: ``'minmax'``（默认）或 ``'std'``，同时决定色带与误差棒的含义。

        Returns:
            汇总输出目录。

        Raises:
            ValueError: 未给 ldir_run，或 band 不是 'minmax'/'std'。
        """
        ldir_run, dir_post_root = self._summary_root(ldir_run, band)
        dir_summary = Path(dir_summary) if dir_summary else dir_post_root
        os.makedirs(dir_summary, exist_ok=True)
        band_label = _BAND_LABEL[band]

        # 1) 学习曲线 -> p_post_rmse_summary（2 panel：上 energy RMSE、下 force RMSE）
        ldf, lrid = self._collect_run_tables(ldir_run, dir_post_root,
                                             Path('training') / 'p_post_learning_curve.txt')
        if ldf:
            df_mean, df_lo, df_hi, df_out = self._ensemble_stats(ldf, _COLS_RMSE, 'epoch', band)
            self._write_table(dir_summary / 'p_post_rmse_summary.txt',
                              ['# Ensemble summary of training-set RMSE vs epoch (from p_post_learning_curve.txt)',
                               '# E_RMSE (meV/atom), F_RMSE (meV/A); physical units',
                               f'# {len(ldf)} run(s): {", ".join(lrid)}',
                               '# n_runs = runs contributing to this epoch; per quantity: _mean/_min/_max/_std '
                               '(std is the sample std, ddof=1, over runs)'],
                              df_out, float_format='16.6f')
            my_plot_epoch_rmse(df_mean, dir_summary / 'p_post_rmse_summary.pdf',
                               df_lo=df_lo, df_hi=df_hi, band_label=band_label.format(n=len(ldf)))
            print(f'  rmse: {len(ldf)} run(s), {len(df_out)} epoch(s) -> p_post_rmse_summary.txt/pdf')
        else:
            print('No run has p_post_learning_curve.txt; skip p_post_rmse_summary.')

        # 2) 按 tag 的 RMSE -> p_post_rmse_by_tag_summary（柱高 = 均值，误差棒 = 上下限）
        ldf, lrid = self._collect_run_tables(ldir_run, dir_post_root,
                                             Path('training') / 'p_post_rmse_by_tag.txt')
        if not ldf:
            print('No run has p_post_rmse_by_tag.txt; skip p_post_rmse_by_tag_summary.')
            return dir_summary

        df_mean, df_lo, df_hi, df_out = self._ensemble_stats(ldf, _COLS_RMSE_BY_TAG, 'tag', band)
        # n_struct / n_fcomp 是数据集属性（各 run 相同），不做统计，按 tag 透传第一个 run 的值
        counts = ldf[0].set_index('tag')[_COLS_RMSE_BY_TAG_COUNT]
        for other, rid in zip(ldf[1:], lrid[1:]):
            if not counts.equals(other.set_index('tag')[_COLS_RMSE_BY_TAG_COUNT]):
                print(f'⚠️ run {rid}: n_struct/n_fcomp differ from run {lrid[0]}; '
                      f'the summary reports run {lrid[0]} counts.')
        for j, c in enumerate(_COLS_RMSE_BY_TAG_COUNT):    # 插在 n_runs 之后、各统计列之前
            df_out.insert(2 + j, c, counts.reindex(df_out['tag'])[c].to_numpy())

        self._write_table(dir_summary / 'p_post_rmse_by_tag_summary.txt',
                          ['# Ensemble summary of per-tag DFT-vs-NNP error (from p_post_rmse_by_tag.txt)',
                           '# E per structure (per atom), F per component; ME = mean(NNP - DFT)',
                           f'# {len(ldf)} run(s): {", ".join(lrid)}',
                           '# n_runs = runs contributing to this tag; n_struct/n_fcomp taken from the first run; '
                           'per quantity: _mean/_min/_max/_std (std is the sample std, ddof=1, over runs)'],
                          df_out, float_format='14.6f')
        my_plot_rmse_by_tag(df_mean, dir_summary / 'p_post_rmse_by_tag_summary.pdf',
                            df_lo=df_lo, df_hi=df_hi, band_label=band_label.format(n=len(ldf)))
        print(f'  rmse_by_tag: {len(ldf)} run(s), {len(df_out)} tag(s) -> p_post_rmse_by_tag_summary.txt/pdf')

        print(f"✅ Ensemble training summary ({band}) -> {dir_summary}")
        return dir_summary


    def post_check_interface_summary(self, ldir_run: list = None, dir_summary: Path = None) -> Path:
        """跨 run 汇总 ``check_interface`` 的一致性判定，写在 ``y_post/`` 下。

        读各 run 的 ``y_post/<run-id>/check_interface/p_post_check_interface.txt``，
        写 ``p_post_check_interface_summary.txt/pdf``：两个 panel（上 energy、下 force），
        横轴为 run，纵轴不写 label，只有两个刻度 —— 1 处标 ``PASS``、2 处标 ``NO``。

        Note:
            纯读取 + 汇总，不触碰任何 run 自己的 ``y_post/<run-id>/`` 产物，可安全重跑。
            某个 run 若因解析失败没有写出 verdict 行，则该点留空（NaN），不会被当成通过。

        Args:
            ldir_run: 参与汇总的 n2p2 训练运行目录列表。
            dir_summary: 汇总输出目录；默认就是 ``y_post/`` 本身。

        Returns:
            汇总输出目录。

        Raises:
            ValueError: 未给 ldir_run。
        """
        ldir_run, dir_post_root = self._summary_root(ldir_run)
        dir_summary = Path(dir_summary) if dir_summary else dir_post_root
        os.makedirs(dir_summary, exist_ok=True)

        rows = []
        for dir_run in ldir_run:
            p = dir_post_root / dir_run.name / 'check_interface' / 'p_post_check_interface.txt'
            if not p.is_file():
                print(f'skip run {dir_run.name}: missing {p}')
                continue
            rows.append({'run': dir_run.name, **self._read_check_interface(p)})
        if not rows:
            print('No run has p_post_check_interface.txt; skip check-interface summary.')
            return dir_summary

        df = pd.DataFrame(rows)[['run', 'epoch', 'tag', 'verdict_E', 'verdict_F',
                                 'delta_E_eV', 'max_dF_comp_eV_A']]
        npass = {c: int((df[c] == VERDICT_CODE['PASS']).sum()) for c in ('verdict_E', 'verdict_F')}
        nfail = {c: int((df[c] == VERDICT_CODE['FAIL']).sum()) for c in ('verdict_E', 'verdict_F')}

        # 表里写人读的 PASS/FAIL，图里用 1/2 编码（纵轴刻度写 PASS/NO）
        df_txt = df.copy()
        code_to_str = {v: k for k, v in VERDICT_CODE.items()}
        for c in ('verdict_E', 'verdict_F'):
            df_txt[c] = df_txt[c].map(lambda v: code_to_str.get(v, 'MISSING'))
        self._write_table(dir_summary / 'p_post_check_interface_summary.txt',
                          ['# Ensemble summary of the nnp-predict (2.3.0) vs LAMMPS hdnnp (2.2.0) interface check',
                           '# verdict PASS/FAIL per run (MISSING = the run wrote no verdict line)',
                           f'# energy {npass["verdict_E"]}/{len(df)} PASS, force {npass["verdict_F"]}/{len(df)} PASS',
                           '# delta_E (eV) targets ~1e-5; max|dF_comp| (eV/A) targets <~1e-5'],
                          df_txt, float_format='16.6e')
        my_plot_check_interface(df, dir_summary / 'p_post_check_interface_summary.pdf')

        verdict = '✅' if not any(nfail.values()) else '❌'
        print(f"{verdict} Ensemble check-interface summary -> {dir_summary}: "
              f"{len(df)} run(s), energy {npass['verdict_E']} PASS / {nfail['verdict_E']} FAIL, "
              f"force {npass['verdict_F']} PASS / {nfail['verdict_F']} FAIL "
              f"(p_post_check_interface_summary.txt/pdf)")
        return dir_summary
