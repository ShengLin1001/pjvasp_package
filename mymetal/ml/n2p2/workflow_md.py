"""
mymetal.ml.n2p2.workflow_md

MD 主动学习（chap_4_TL.pdf Stage I）的调度器：把「已训练好的一组势函数」变成
「一批值得送 DFT 标记的新构型」，再回灌训练集开始下一轮训练。

``PeiN2p2MD`` 继承 :class:`mymetal.ml.n2p2.workflow.PeiN2p2`，只新增 MD 主动学习这条
支线；数据/对称函数/训练/后处理仍走父类原有方法（``submit_train`` / ``post_training`` /
``post_properties`` / ``post_epoch_scan`` …），不做任何重写。父类已经很长，MD 相关的
目录约定、LAMMPS 采样、committee 一致性筛选与 DFT 标记全部落在本模块。

一轮（iround）的完整链路：

    collect_potentials      train/.../y_n2p2_train/y_dir/<run>/ 的最后 n_last 条势
                            -> potentials/<run>/<epoch>/
    submit_md               每条势一个 LAMMPS 升温 MD 作业
                            -> data_md/<iround>/y_dir/<run>_<epoch>/
    collect_md_frames       所有轨迹的帧汇成一份 n2p2 数据集
                            -> data_md/<iround>/y_frames/input.data.frames
    submit_compare          每条势对**整条**帧集跑一次 nnp-dataset
                            -> data_md/<iround>/y_compare/y_dir/<run>_<epoch>/energy.comp.*
    collect_compare         逐帧汇总各势预测 -> y_compare/y_compare_data_md.txt
    select_md_structures    按 committee 发散度挑结构 -> y_compare/input.data.md
    prepare_dft             逐结构 VASP 静态 -> data_dft/<iround>/y_dir/<i>/
    collect_dft             读回 OUTCAR -> data_dft/<iround>/y_dir/input.data.dft
    append_to_train         input.data.ini + 各轮 input.data.dft -> data/train/input.data
    （随后照常 PeiN2p2.submit_train(dir_run=train/<iround>/y_n2p2_train/y_dir/<run>)）

为什么用 ``nnp-dataset`` 而不是 ``nnp-predict``：nnp-predict 一次只算 input.data 里的
**一个**结构，逐帧逐势建目录会变成「帧数 × 势数」个工作目录（1 万帧 × 120 条势 = 120 万个
目录、上百 GB 势函数副本），在共享文件系统上不可行。nnp-dataset 一次吃完整份 input.data，
按结构逐行写 ``energy.comp.%04d``（index / N / Eref_phys / Ennp_phys），信息完全等价，
目录数降到「势数」量级。

目录约定刻意贴着 ``pei_slurm_univ_submit`` 的发现规则：该引擎递归找出任意深度的每个
``y_dir``，把它的**一级**子目录当作作业目录。所以 MD 与 compare 的作业目录都平铺成
``y_dir/<run>_<epoch>/``，而不是再嵌一层。

Classes:
    PeiN2p2MD: n2p2 MD 主动学习（采样 -> committee 筛选 -> DFT 标记 -> 回灌）调度器。

Functions:
    read_lammps_dump: 读 LAMMPS custom dump（含三斜盒子换算、每原子力与每原子能量）。
    read_temp_ladder: 读 md_heating.mod 写出的温度阶梯表。
"""

from mymetal.ml.n2p2.workflow import PeiN2p2
from mymetal.ml.n2p2.dataset import nnpdata
from mymetal.ml.n2p2.calculate.sf import get_largest_rc_from_input_nn

from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
import numpy as np
import pandas as pd
import os
from pathlib import Path
import shutil


# ===== MD 采样的默认参数 =====
# 温度阶梯默认沿用 chap_4_TL.pdf Stage I：25 -> 450 K、25 K 一档、共 18 档。
# ⚠️ 该设置来自 Mg（Tm≈923 K，450 K ≈ 0.49 Tm）；Au 的 Tm≈1337 K，同比例应到 ~650 K。
#    CLAUDE.md 明确要求 Mg 示例只迁移方法框架、不照搬数值，故这里把它做成可调参数，
#    默认值留在 LT 口径上便于对照，实际 Au 训练时应按上面的比例重标定并记录在 README。
#
# 超胞：FCC 常规胞 4 原子 -> 2x2x2 = 32；HCP 六方胞 2 原子 -> 3x3x2 = 36。
# 这些帧最终要逐个送 DFT 静态标记，原子数直接决定 DFT 代价，也要与训练集口径相称
# （stage1/0 把训练结构规整到 [20,30] 逼近 24）。
#
# 采样节奏：每档 n_equil 步平衡（不采样）+ n_prod 步采样、每 n_dump 步一帧
# -> 每档 20 帧，单相 360 帧，单个作业（FCC+HCP）720 帧。
DICT_MD_PARAMS_DEFAULT = {
    'lphase': ['fcc', 'hcp'],
    'llatid': [2, 1],                    # stretch_model.mod 的 lat 编号：1=hcp 2=fcc 3=bcc
    'laa0':   [4.08, 2.88],              # 初始晶格常数猜测（Å），升温前先 box/relax
    'lnx':    [2, 3],
    'lny':    [2, 3],
    'lnz':    [2, 2],
    't_start': 25.0,
    't_step':  25.0,
    'n_temp':  18,
    'n_equil': 2000,
    'n_prod':  10000,
    'n_dump':  500,
    'dt':      0.001,                    # metal 单位为 ps
    't_damp':  0.1,
    'p_damp':  1.0,
}

# md_params.mod 里 index 风格（逐相一个取值）与 equal 风格（全局标量）的键分组。
# 两组必须分开渲染：index 变量靠 next 推进相循环，equal 变量每相重新求值。
_LKEY_MD_INDEX = ['lphase', 'llatid', 'laa0', 'lnx', 'lny', 'lnz']
_LKEY_MD_SCALAR = ['t_start', 't_step', 'n_temp', 'n_equil', 'n_prod', 'n_dump',
                   'dt', 't_damp', 'p_damp']
# index 变量在 md_params.mod / md_heating_template.in 里的 LAMMPS 变量名
_DICT_MD_INDEX_NAME = {'lphase': 'phase', 'llatid': 'latid', 'laa0': 'aa0',
                       'lnx': 'nx', 'lny': 'ny', 'lnz': 'nz'}

# 各产物表的列顺序（写表与读回共用一份定义，避免两边漂移）
_COLS_FRAMES = ['index', 'run', 'epoch', 'phase', 'T_set_K', 'step', 'natoms', 'E_md_eV']
_COLS_COMPARE_STAT = ['index', 'phase', 'T_set_K', 'natoms',
                      'E_mean_eV/at', 'E_std_meV/at', 'E_spread_meV/at', 'n_pot']
_COLS_SELECT = ['rank', 'index', 'phase', 'T_set_K', 'natoms',
                'E_std_meV/at', 'E_spread_meV/at', 'tag']

# LAMMPS dump 的列名 -> 本模块内部字段名。改 md_heating.mod 的 dump 列时这里同步改。
_DUMP_COLS_REQUIRED = ['id', 'type', 'x', 'y', 'z', 'fx', 'fy', 'fz', 'c_peatom1']


def read_temp_ladder(path_temp: Path = None) -> list:
    """读 ``md_heating.mod`` 写出的温度阶梯表 ``y_md_<phase>_temp.txt``。

    表头是 ``k T_set step_begin step_end``，其后每档一行。step 区间用来给帧标温度：
    dump 与温度记录不是同一时刻落盘的，逐帧精确匹配会错位，区间归属则不会。

    Args:
        path_temp: y_md_<phase>_temp.txt 路径。

    Returns:
        ``[(T_set, step_begin, step_end), ...]``，按 step_begin 升序。

    Raises:
        FileNotFoundError: 文件不存在或为空。
    """
    path_temp = Path(path_temp)
    if not (path_temp.is_file() and path_temp.stat().st_size > 0):
        raise FileNotFoundError(f"❌ Missing/empty temperature ladder: {path_temp}")

    lladder = []
    for line in path_temp.read_text(encoding='utf-8').splitlines():
        lfield = line.split()
        if len(lfield) != 4 or lfield[0] == 'k':      # 表头行
            continue
        lladder.append((float(lfield[1]), int(lfield[2]), int(lfield[3])))
    if not lladder:
        raise FileNotFoundError(f"❌ No temperature rows parsed from {path_temp}")
    return sorted(lladder, key=lambda row: row[1])


def _bounds_to_cell(lbound: list, ltilt: list) -> np.ndarray:
    """把 LAMMPS dump 的 BOX BOUNDS（bound 值 + 倾斜量）还原成 3x3 晶胞。

    dump 写的是把倾斜量算进去之后的**外接**边界，必须按 LAMMPS 文档反解回 lo/hi
    才能得到真实盒长；直接拿 bound 相减会让三斜胞（HCP 的 xy 倾斜）整体偏大。

    Args:
        lbound: ``[(xlo_b, xhi_b), (ylo_b, yhi_b), (zlo_b, zhi_b)]``。
        ltilt: ``[xy, xz, yz]``；正交胞传 ``[0.0, 0.0, 0.0]``。

    Returns:
        行向量形式的 3x3 晶胞矩阵。
    """
    (xlo_b, xhi_b), (ylo_b, yhi_b), (zlo_b, zhi_b) = lbound
    xy, xz, yz = ltilt

    xlo = xlo_b - min(0.0, xy, xz, xy + xz)
    xhi = xhi_b - max(0.0, xy, xz, xy + xz)
    ylo = ylo_b - min(0.0, yz)
    yhi = yhi_b - max(0.0, yz)

    return np.array([[xhi - xlo, 0.0, 0.0],
                     [xy, yhi - ylo, 0.0],
                     [xz, yz, zhi_b - zlo_b]], dtype=float)


def read_lammps_dump(path_dump: Path = None, lele: list = None) -> list:
    """读 ``md_heating.mod`` 写出的 custom dump，逐帧还原成带能量/力的 ASE Atoms。

    自己解析而不是走 ``ase.io.read(format='lammps-dump-text')``：本 dump 多带一列
    每原子能量 ``c_peatom1``，ASE 的读取器不认识这个列名会直接丢掉，而总能量正是
    committee 对比与 input.data 落盘所必需的（逐帧总能 = 该列求和）。

    Args:
        path_dump: y_md_<phase>.lammpstrj 路径。
        lele: LAMMPS 原子类型序（type i -> lele[i-1]），一般取 ``self.lele``。

    Returns:
        ``[(timestep, Atoms), ...]``，Atoms 上挂好 energy/forces 的 SinglePointCalculator。

    Raises:
        FileNotFoundError: 文件不存在或为空。
        ValueError: dump 缺少必需列，或原子数与实际行数不一致。
    """
    path_dump = Path(path_dump)
    if not (path_dump.is_file() and path_dump.stat().st_size > 0):
        raise FileNotFoundError(f"❌ Missing/empty dump: {path_dump}")
    if not lele:
        raise ValueError("❌ lele must be given to map LAMMPS types back to elements.")

    lline = path_dump.read_text(encoding='utf-8').splitlines()
    lframe = []
    i, n = 0, len(lline)
    while i < n:
        if not lline[i].startswith('ITEM: TIMESTEP'):
            i += 1
            continue

        timestep = int(lline[i + 1].split()[0])
        natom = int(lline[i + 3].split()[0])

        # "ITEM: BOX BOUNDS xy xz yz pp pp pp"：带 xy/xz/yz 才是三斜胞，每行多一个倾斜量
        head_box = lline[i + 4]
        if_tri = 'xy' in head_box
        lbound, ltilt = [], []
        for j in range(3):
            lfield = [float(x) for x in lline[i + 5 + j].split()]
            lbound.append((lfield[0], lfield[1]))
            ltilt.append(lfield[2] if if_tri else 0.0)
        cell = _bounds_to_cell(lbound, ltilt)

        # "ITEM: ATOMS id type x y z fx fy fz c_peatom1"
        i_atoms = i + 8
        lcol = lline[i_atoms].split()[2:]
        lmiss = [c for c in _DUMP_COLS_REQUIRED if c not in lcol]
        if lmiss:
            raise ValueError(f"❌ {path_dump} frame at step {timestep} misses dump column(s) "
                             f"{lmiss}; got {lcol}. Check md_heating.mod's dump line.")
        dict_icol = {name: k for k, name in enumerate(lcol)}

        arr = np.array([lline[i_atoms + 1 + k].split() for k in range(natom)], dtype=float)
        if arr.shape[0] != natom:
            raise ValueError(f"❌ {path_dump} frame at step {timestep}: expected {natom} atom "
                             f"rows, parsed {arr.shape[0]}.")

        ltype = arr[:, dict_icol['type']].astype(int)
        if ltype.max() > len(lele):
            raise ValueError(f"❌ {path_dump}: atom type {ltype.max()} exceeds element list {lele}.")
        positions = arr[:, [dict_icol['x'], dict_icol['y'], dict_icol['z']]]
        forces = arr[:, [dict_icol['fx'], dict_icol['fy'], dict_icol['fz']]]
        energy = float(arr[:, dict_icol['c_peatom1']].sum())

        atoms = Atoms(symbols=[lele[t - 1] for t in ltype],
                      positions=positions, cell=cell, pbc=True)
        atoms.calc = SinglePointCalculator(atoms, energy=energy, forces=forces)
        lframe.append((timestep, atoms))

        i = i_atoms + 1 + natom
    return lframe


class PeiN2p2MD(PeiN2p2):
    """n2p2 MD 主动学习调度器（chap_4_TL.pdf Stage I 的 Au 版实现）。

    继承 :class:`~mymetal.ml.n2p2.workflow.PeiN2p2` 的全部目录约定、Slurm 提交/等待
    工具（``_sbatch_capture`` / ``_snapshot_jobids`` / ``wait_jobs``）、重建守卫
    ``_guard_rebuild_dir`` 与完成阀门 ``_skip_done_stage``，只新增 MD 采样 ->
    committee 筛选 -> DFT 标记 -> 回灌训练集这一条支线。训练本身仍调父类
    ``submit_train(dir_run=train/<iround>/y_n2p2_train/y_dir/<run>)``，不重写。

    Attributes:
        dir_potentials: 从上一轮训练收集来的势函数库根目录。
        dir_data_md: 逐轮 MD 采样与 committee 对比的根目录。
        dir_data_dft: 逐轮 DFT 标记的根目录。

    TODO:
        - 目前建模只支持单元素（stretch_model.mod 的 ``create_box 1 box``）。
          多元素体系要先扩展该模板，本类的 pair_coeff / mass 侧已按 self.lele 泛化。
    """

    def __init__(self, dir_potentials: Path = Path('./potentials'),
                 dir_data_md: Path = Path('./data_md'),
                 dir_data_dft: Path = Path('./data_dft'),
                 **kwargs):
        """初始化 MD 主动学习目录，其余参数原样透传给父类。

        Args:
            dir_potentials: 势函数库根目录（``<run>/<epoch>/`` 两级）。
            dir_data_md: MD 采样与 committee 对比根目录（一级子目录 = 轮次）。
            dir_data_dft: DFT 标记根目录（一级子目录 = 轮次）。
            **kwargs: 透传给 :class:`PeiN2p2`（``dir_root`` / ``dir_lmp_utils`` /
                ``path_python`` 等必填项仍由调用方给出）。
        """
        super().__init__(**kwargs)

        self.dir_potentials = Path(dir_potentials)
        self.dir_data_md = Path(dir_data_md)
        self.dir_data_dft = Path(dir_data_dft)

        os.makedirs(self.dir_potentials, exist_ok=True)
        os.makedirs(self.dir_data_md, exist_ok=True)
        os.makedirs(self.dir_data_dft, exist_ok=True)


    # ===== 阶段 i：收集势函数 =====
    @staticmethod
    def _resolve_lepoch(dir_run: Path = None, lepoch: list = None,
                        n_last: int = None) -> list:
        """定出单个训练 run 要取哪几个 epoch 的权重快照。

        两种选法，``lepoch`` 优先：

        - ``lepoch=[1000]`` / ``lepoch=[500, 600, 700, 800]``：**按 epoch 号显式点名**。
          这是默认用法 —— committee 的独立性来自不同的 ``random_seed``，即「不同的 run」；
          同一个 run 的相邻 epoch 是同一条优化轨迹上的快照，彼此高度相关，不是独立样本
          （实测 run 内 std 只有 run 间 std 的 1/7）。想扩大 committee 就加 run，
          想沿训练轨迹看势函数怎么演化才多点几个 epoch。
        - ``n_last=20``：取末尾若干个（旧行为，保留给「扫末段收敛性」这类用途）。

        Args:
            dir_run: 训练运行目录（含 ``weights.<Z>.<epoch>.out``）。
            lepoch: 要取的 epoch 号列表，接受 int（``1000``）或 str（``'001000'``）；
                按 epoch **数值**匹配，不依赖零填充宽度。
            n_last: 仅当 ``lepoch`` 为 None 时生效；取末尾多少个，None 或 <=0 表示全部。

        Returns:
            epoch 字符串列表（权重文件里的原始零填充形式），按 epoch 升序。

        Raises:
            FileNotFoundError: 目录里没有任何权重文件。
            ValueError: 点名的 epoch 在该 run 里不存在（会列出可用范围，避免
                「少拷了一条势却一路跑到 compare 才发现列数对不上」）。
        """
        dir_run = Path(dir_run)
        lepoch_avail = sorted({w.name.split('.')[-2] for w in dir_run.glob('weights.*.*.out')})
        if not lepoch_avail:
            raise FileNotFoundError(f"❌ No weights.*.*.out in {dir_run}. Training not finished?")

        if lepoch is not None:
            # 按数值建索引：训练里 epoch 段是 6 位零填充，但调用方写 1000 更自然。
            dict_epoch = {int(e): e for e in lepoch_avail}
            lout, lmiss = [], []
            for e in lepoch:
                key = int(e)
                (lout.append(dict_epoch[key]) if key in dict_epoch else lmiss.append(key))
            if lmiss:
                raise ValueError(
                    f"❌ epoch(s) {lmiss} have no weights.*.<epoch>.out in {dir_run}. "
                    f"Available: {int(lepoch_avail[0])} .. {int(lepoch_avail[-1])} "
                    f"({len(lepoch_avail)} snapshot(s)). Check `write_weights_epoch` "
                    f"in input.nn — weights are only dumped every N epochs.")
            return sorted(set(lout), key=int)

        if n_last and n_last > 0:
            return lepoch_avail[-n_last:]
        return lepoch_avail

    def collect_potentials(self, dir_train_src: Path = None, lepoch: list = None,
                           n_last: int = None, lexclude: list = None,
                           if_force_rebuild: bool = False,
                           if_skip_if_done: bool = True) -> list:
        """把上一轮训练的若干条势函数收进 ``potentials/<run>/<epoch>/``。

        每个 ``<epoch>/`` 是一份**自包含**的势函数三件套（input.nn + scaling.data +
        每元素一份 ``weights.<Z>.data``），nnp-predict / nnp-dataset / LAMMPS hdnnp
        都按 cwd 或 ``dir`` 直接读，无需再回训练目录取文件。

        Note:
            目录名用**源 run 的 basename 与 epoch**（如 ``001/001000``）而不是 0..5 / 0..19
            的顺序编号：``lexclude`` 一旦排掉某个 run，顺序编号就与训练目录对不上，而
            committee 的每一列都要能追回是哪个 run 的哪个 epoch。

        Args:
            dir_train_src: 上一轮的 ``y_n2p2_train`` 目录（其下为 ``y_dir/<run>/``）。
            lepoch: 每个 run 取哪几个 epoch，如 ``[1000]``（默认口径，一 run 一条）或
                ``[500, 600, 700, 800]``。接受 int 或零填充 str，按数值匹配。
                与 ``n_last`` 二选一，本参数优先；两者都不给时默认 ``[1000]``。
            n_last: 旧口径 —— 每个 run 取末尾多少条。仅当 ``lepoch`` 为 None 时生效。
                ⚠️ 用它构造 committee 是无效的：同 run 相邻 epoch 高度相关，
                样本量涨了但有效样本量还是 run 数（见 :meth:`_resolve_lepoch`）。
            lexclude: 要排除的 run basename 列表（如训练种子不好的那几个）。
            if_force_rebuild: ``potentials/`` 非空时是否强制清空重建。
            if_skip_if_done: 完成阀门；清单文件已非空则直接读回、不重拷。

        Returns:
            ``[(run, epoch), ...]``，按 (run, epoch) 升序。

        Raises:
            FileNotFoundError: 找不到 ``y_dir``、没有可用 run，或某个 run 缺权重/势文件。
        """
        os.chdir(self.dir_root)
        if dir_train_src is None:
            raise ValueError("❌ dir_train_src must be given (the previous y_n2p2_train dir).")
        if lepoch is None and not n_last:
            lepoch = [1000]          # 默认：每个 run 一条末期势，committee 大小 = run 数
        dir_train_src = Path(dir_train_src)
        path_manifest = self.dir_potentials / 'p_post_potentials.txt'

        if self._skip_done_stage([path_manifest], f'collect_potentials {self.dir_potentials}',
                                 if_skip_if_done=if_skip_if_done):
            return self._read_lpotential()

        dir_ydir_src = dir_train_src / 'y_dir'
        if not dir_ydir_src.is_dir():
            raise FileNotFoundError(f"❌ Missing {dir_ydir_src}.")

        lexclude = set(lexclude or [])
        ldir_run = sorted(d for d in dir_ydir_src.iterdir()
                          if d.is_dir() and d.name not in lexclude)
        if not ldir_run:
            raise FileNotFoundError(f"❌ No usable run dir under {dir_ydir_src} "
                                    f"(excluded: {sorted(lexclude)}).")
        if lexclude:
            print(f"  ⚠️  excluding {len(lexclude)} run dir(s): {sorted(lexclude)}")

        self._guard_rebuild_dir(self.dir_potentials, if_force_rebuild=if_force_rebuild,
                                label='potentials')
        os.makedirs(self.dir_potentials, exist_ok=True)

        lpotential = []
        for dir_run in ldir_run:
            for fn in ['input.nn', 'scaling.data']:
                if not self._is_file_nonempty(dir_run / fn):
                    raise FileNotFoundError(f"❌ Missing/empty {dir_run / fn}.")
            lepoch_run = self._resolve_lepoch(dir_run, lepoch=lepoch, n_last=n_last)
            for epoch in lepoch_run:
                dir_pot = self.dir_potentials / dir_run.name / epoch
                os.makedirs(dir_pot, exist_ok=True)
                shutil.copy(dir_run / 'input.nn', dir_pot / 'input.nn')
                shutil.copy(dir_run / 'scaling.data', dir_pot / 'scaling.data')
                lw = sorted(dir_run.glob(f'weights.*.{epoch}.out'))
                if not lw:
                    raise FileNotFoundError(f"❌ No weights.*.{epoch}.out in {dir_run}.")
                for w in lw:
                    # hdnnp / nnp-* 一律按 weights.<Z>.data 找权重，去掉中间的 epoch 段
                    shutil.copy(w, dir_pot / f"weights.{w.name.split('.')[1]}.data")
                lpotential.append((dir_run.name, epoch))
            print(f"  ✅ {dir_run.name}: {len(lepoch_run)} potential(s) "
                  f"({lepoch_run[0]} .. {lepoch_run[-1]})")

        df = pd.DataFrame(lpotential, columns=['run', 'epoch'])
        self._write_table(path_manifest,
                          ['# potentials collected for MD active learning',
                           f'# source     {dir_train_src}',
                           f'# lepoch     {lepoch if lepoch is not None else "(last %s)" % n_last}',
                           f'# committee  {len(ldir_run)} run(s); the independent dimension '
                           f'is the run (random_seed), not the epoch',
                           f'# excluded   {sorted(lexclude) if lexclude else "(none)"}',
                           f'# layout     potentials/<run>/<epoch>/'
                           f'{{input.nn, scaling.data, weights.<Z>.data}}'],
                          df, float_format='16.8f')
        print(f"📊 collect_potentials: {len(lpotential)} potential(s) from {len(ldir_run)} run(s) "
              f"-> {self.dir_potentials}")
        return lpotential

    def _read_lpotential(self) -> list:
        """读回 ``potentials/p_post_potentials.txt`` 的 (run, epoch) 清单。

        Returns:
            ``[(run, epoch), ...]``。

        Raises:
            FileNotFoundError: 清单文件不存在或为空（未跑过 collect_potentials）。
        """
        path_manifest = self.dir_potentials / 'p_post_potentials.txt'
        if not self._is_file_nonempty(path_manifest):
            raise FileNotFoundError(f"❌ Missing {path_manifest}. Run collect_potentials first.")
        df = pd.read_csv(path_manifest, sep=r'\s+', comment='#', dtype=str)
        return [(row['run'], row['epoch']) for _, row in df.iterrows()]

    @staticmethod
    def _potential_id(run: str, epoch: str) -> str:
        """势函数在作业目录名与表格列名里的统一标识（``<run>_<epoch>``）。

        Args:
            run: 训练 run 的 basename。
            epoch: 权重快照的 epoch 字符串。

        Returns:
            标识字符串。
        """
        return f'{run}_{epoch}'


    # ===== 三个提交入口共用的部件 =====
    @staticmethod
    def _clamp_chunks(dict_args: dict = None, n_job: int = None) -> dict:
        """把 ``chunks`` 夹到本次实际提交的作业目录数以内。

        ``pei_slurm_univ_submit`` 明确拒绝 ``chunks > 作业目录数``（直接 fail 退出）。
        首轮一般不会踩到（作业上百个），但**补投**时只剩几个未完成目录，仍按 chunks=10
        提交就会整批失败——而这恰恰是最需要它跑起来的时候。

        Args:
            dict_args: 传给提交器的参数字典。
            n_job: 本次真正要提交的作业目录数。

        Returns:
            夹取后的新字典（不改原字典）。
        """
        dict_args = dict(dict_args or {})
        n_job = max(1, int(n_job))
        chunks = int(dict_args.get('chunks', 1))
        if chunks > n_job:
            print(f"  ⚠️  chunks {chunks} > {n_job} job dir(s) to submit; clamped to {n_job} "
                  f"(the submit engine rejects chunks > job count).")
            dict_args['chunks'] = n_job
        return dict_args

    # ===== 阶段 ii：逐条势函数做升温 MD =====
    @staticmethod
    def _render_md_params(dict_md_params: dict = None) -> str:
        """把 MD 参数字典渲染成 ``md_params.mod`` 的内容。

        index 风格变量（逐相一个取值）与 equal 风格标量分开渲染：前者靠
        ``md_heating_template.in`` 的 ``next`` 推进相循环，长度必须彼此相等；
        后者每相重新求值。用普通字符串拼接而不是 f-string —— 生成的是 LAMMPS
        脚本，``${...}`` 与 f-string 的花括号会互相打架。

        Args:
            dict_md_params: 见 :data:`DICT_MD_PARAMS_DEFAULT`。

        Returns:
            md_params.mod 的完整文本。

        Raises:
            ValueError: 缺键，或各 index 列表长度不一致。
        """
        dict_md_params = {**DICT_MD_PARAMS_DEFAULT, **(dict_md_params or {})}

        lmiss = [k for k in _LKEY_MD_INDEX + _LKEY_MD_SCALAR if k not in dict_md_params]
        if lmiss:
            raise ValueError(f"❌ dict_md_params misses key(s): {lmiss}")
        lnphase = {k: len(dict_md_params[k]) for k in _LKEY_MD_INDEX}
        if len(set(lnphase.values())) != 1:
            raise ValueError(f"❌ per-phase lists must have equal length, got {lnphase}. "
                             f"md_heating_template.in advances them with one `next`.")

        lline = ['# md_params.mod — 自动生成，请勿手改（每次 submit_md 覆盖）',
                 '# 相列表用 index 风格变量：LAMMPS 的 clear 不清除 input script variables，',
                 '# 故 md_heating_template.in 能用一个 in 文件跨 clear 依次跑完所有相。',
                 '']
        for key in _LKEY_MD_INDEX:
            lvalue = ' '.join(str(v) for v in dict_md_params[key])
            lline.append('variable ' + _DICT_MD_INDEX_NAME[key] + ' index ' + lvalue)
        lline.append('')
        lline.append('# 温度阶梯 T_k = t_start + (k-1)*t_step，共 n_temp 档；逐档升温、不重置速度')
        for key in _LKEY_MD_SCALAR:
            lline.append('variable ' + key + ' equal ' + str(dict_md_params[key]))
        lline.append('')
        return '\n'.join(lline)

    def _md_job_markers(self, dir_job: Path = None, lphase: list = None) -> list:
        """单个 MD 作业目录里 :meth:`collect_md_frames` 会去读的全部产物。

        Args:
            dir_job: ``data_md/<iround>/y_dir/<run>_<epoch>/``。
            lphase: 相名列表。

        Returns:
            产物文件路径列表。
        """
        dir_job = Path(dir_job)
        lpath = []
        for phase in lphase:
            lpath.append(dir_job / f'y_md_{phase}.lammpstrj')
            lpath.append(dir_job / f'y_md_{phase}_temp.txt')
        return lpath

    def _prepare_md_job(self, dir_job: Path = None, run: str = None, epoch: str = None,
                        md_params_text: str = None) -> None:
        """为单条势函数铺好一个 MD 作业目录（仅脚手架，不做 sed）。

        与 :meth:`PeiN2p2._prepare_epoch_dir` 同款分工：Python 只铺 potential/ +
        template/ + md_params.mod，所有 sed 与运行交给作业里的
        ``pei_lmp_run_md_heating``。

        Args:
            dir_job: 目标作业目录。
            run: 势函数来源的训练 run basename。
            epoch: 势函数的 epoch 字符串。
            md_params_text: 渲染好的 md_params.mod 内容。

        Raises:
            FileNotFoundError: 势函数目录缺文件，或 lmp_utils 模板目录不存在。
        """
        dir_job = Path(dir_job)
        dir_pot_src = self.dir_potentials / run / epoch
        for fn in ['input.nn', 'scaling.data']:
            if not self._is_file_nonempty(dir_pot_src / fn):
                raise FileNotFoundError(f"❌ Missing/empty {dir_pot_src / fn}. "
                                        f"Run collect_potentials first.")
        lw = sorted(dir_pot_src.glob('weights.*.data'))
        if not lw:
            raise FileNotFoundError(f"❌ No weights.*.data in {dir_pot_src}.")

        dir_potential = dir_job / 'potential'
        os.makedirs(dir_potential, exist_ok=True)
        for path_src in [dir_pot_src / 'input.nn', dir_pot_src / 'scaling.data'] + lw:
            shutil.copy(path_src, dir_potential / path_src.name)

        # 整目录拷一份模板到作业里：runner 就地 sed pair_style/pair_coeff 与
        # stretch_model.mod 的 lat/aa 占位符，绝不回头改 lmp_utils 里的源文件。
        shutil.copytree(self.dir_lmp_utils / 'template', dir_job / 'template', dirs_exist_ok=True)
        (dir_job / 'md_params.mod').write_text(md_params_text, encoding='utf-8')

    def submit_md(self, iround: int = 0, lpotential: list = None,
                  dict_md_params: dict = None, if_sbatch: bool = False,
                  dict_args_to_submit: dict = None,
                  if_force_rebuild: bool = False, if_skip_if_done: bool = True) -> Path:
        """准备并可选提交本轮的全部升温 MD 作业（每条势函数一个）。

        作业目录平铺在 ``data_md/<iround>/y_dir/<run>_<epoch>/``：
        ``pei_slurm_univ_submit`` 把每个 ``y_dir`` 的**一级**子目录当作业目录，多嵌一层
        就会把「按 run 分组的中间目录」当成作业本身投出去。

        Note:
            随机种子刻意**不在这里固定**：由作业里的 ``pei_lmp_run_md_heating`` 每次运行
            现取一个写进 md_seed.mod。种子写死会让同一条势函数重跑得到完全相同的轨迹，
            committee 采样就失去了独立性。

        Args:
            iround: MD 扩充轮次（0, 1, 2, ...）。
            lpotential: ``[(run, epoch), ...]``；为 None 时读 potentials 的清单。
            dict_md_params: 覆盖 :data:`DICT_MD_PARAMS_DEFAULT` 的 MD 参数。
            if_sbatch: 是否真正提交 Slurm 作业。
            dict_args_to_submit: 传给通用 Slurm 提交器的参数字典。
            if_force_rebuild: 作业目录非空时是否强制清空重建。
            if_skip_if_done: 完成阀门，逐作业生效：产物齐全的目录不重铺也不重投；
                脚手架已铺但没结果的目录原样保留、照常提交（补投路径）。

        Returns:
            本轮 MD 根目录 ``data_md/<iround>/``。

        Raises:
            FileNotFoundError: 缺势函数、缺 lmp_utils 模板目录。
        """
        os.chdir(self.dir_root)
        if dict_args_to_submit is None:
            # 单核单节点：一条 MD 轨迹的体系只有几十个原子，多核反而被通信吃掉；
            # chunks 把上百个作业压成少量并发调度流，避免一次性冲垮 slurmctld。
            dict_args_to_submit = {'preset': 'zcm6_lammps_0', 'mode': 'each_subdir',
                                   'chunks': 10, 'nodes': 1, 'ncores': 1,
                                   'launcher_type': 'srun', 'if_use_my_launcher': True}
        if not self.dir_lmp_utils.is_dir():
            raise FileNotFoundError(f"❌ Missing lmp_utils source dir {self.dir_lmp_utils}.")

        lpotential = lpotential if lpotential is not None else self._read_lpotential()
        dict_md_params = {**DICT_MD_PARAMS_DEFAULT, **(dict_md_params or {})}
        md_params_text = self._render_md_params(dict_md_params)
        lphase = dict_md_params['lphase']

        dir_round = self.dir_data_md / str(iround)
        dir_ydir = dir_round / 'y_dir'
        os.makedirs(dir_ydir, exist_ok=True)

        ldir_job = [dir_ydir / self._potential_id(run, epoch) for run, epoch in lpotential]
        ldir_done = ([d for d in ldir_job
                      if all(self._is_file_nonempty(p) for p in self._md_job_markers(d, lphase))]
                     if if_skip_if_done else [])
        if if_skip_if_done and len(ldir_done) == len(ldir_job):
            print(f"⏭️ skip submit_md {dir_round}: all {len(ldir_job)} MD job dir(s) already carry "
                  f"complete trajectories. Pass if_skip_if_done=False to redo them.")
            self.last_jobids = []
            return dir_round

        n_kept = 0
        for (run, epoch), dir_job in zip(lpotential, ldir_job):
            if if_skip_if_done and self._is_file_nonempty(dir_job / 'md_params.mod') \
                    and self._dir_covers_source(dir_job / 'template', self.dir_lmp_utils / 'template'):
                n_kept += 1
                continue
            self._guard_rebuild_dir(dir_job, if_force_rebuild=if_force_rebuild,
                                    label=f'md job {dir_job.name}')
            self._prepare_md_job(dir_job, run, epoch, md_params_text)
        print(f"Prepared {len(ldir_job) - n_kept} MD job dir(s) under {dir_ydir} "
              f"(kept {n_kept} already-prepared, {len(ldir_done)} of them with complete results).")

        # runner 的 4 个位置参数（第 4 个 lmp_template_path 用默认 ./template，不在此传）：
        #   $1 pair_style   hdnnp 整串；dir 用运行期 $(pwd)/potential（逐作业不同）
        #   $2 pair_coeff   元素串（hdnnp 用元素名，不是势文件路径）
        #   $3 mass_content 写入 general_mass.mod 的多行内容（\n 由 runner 的 printf '%b' 展开）
        cutoff = self._get_md_cutoff(lpotential)
        from ase.data import atomic_numbers, atomic_masses
        pair_style = (f'hdnnp {cutoff:.2f} dir $(pwd)/potential '
                      'showew no showewsum 100 resetew yes maxew 1000000 cflength 1.0 cfenergy 1.0')
        pair_coeff = ' '.join(self.lele)
        mass_content = '\\n'.join(f'mass {i + 1} {atomic_masses[atomic_numbers[el]]:.4f}'
                                  for i, el in enumerate(self.lele))
        cmd = f'pei_lmp_run_md_heating "{pair_style}" "{pair_coeff}" "{mass_content}"'

        # 只把**未完成**的作业交给引擎：不给 -lsubdir 会扫遍整个 y_dir，把已跑完的轨迹重投覆盖。
        ldir_todo = [d for d in ldir_job if d not in set(ldir_done)]
        args_to_submit = {**self._clamp_chunks(dict_args_to_submit, len(ldir_todo)),
                          'cmd': cmd,
                          'path_root': str(dir_round.resolve()),
                          'dir_root': './y_dir',
                          'if_sbatch': True if if_sbatch else False}
        if ldir_done:
            args_to_submit['lsubdir'] = [d.name for d in ldir_todo]
            print(f"  submit only {len(ldir_todo)} incomplete MD job(s) "
                  f"(skip {len(ldir_done)} already complete)")
        cli_args = self._dict_to_cli_args(args_to_submit)

        self.last_jobids = []
        if if_sbatch:
            # pei_slurm_univ_submit 不回传作业号；用提交前后的 squeue 差集捕获新增作业。
            before = self._snapshot_jobids()
            os.system(f'pei_slurm_univ_submit.py {cli_args}')
            after = self._snapshot_jobids()
            self.last_jobids = sorted(after - before)
            print(f"  submit_md observed {len(self.last_jobids)} new Slurm job(s): {self.last_jobids}")
        return dir_round

    def _get_md_cutoff(self, lpotential: list = None) -> float:
        """pair_style hdnnp 的近邻表 cutoff：最大对称函数 rc 再留 0.01 Å 余量。

        Args:
            lpotential: ``[(run, epoch), ...]``；只读第一条（同一轮的势共用一套 SF）。

        Returns:
            cutoff（Å）。

        Raises:
            FileNotFoundError: 势函数目录里没有 input.nn。
        """
        run, epoch = lpotential[0]
        path_input_nn = self.dir_potentials / run / epoch / 'input.nn'
        if not self._is_file_nonempty(path_input_nn):
            raise FileNotFoundError(f"❌ Missing/empty {path_input_nn}.")
        return get_largest_rc_from_input_nn(str(path_input_nn)) + 0.01


    # ===== 阶段 iii-a：把所有轨迹汇成一份帧集 =====
    @staticmethod
    def _temperature_of_step(step: int = None, lladder: list = None) -> float:
        """按 step 区间给单帧标温度。

        Args:
            step: 帧的 LAMMPS timestep。
            lladder: :func:`read_temp_ladder` 的返回值。

        Returns:
            该帧所属温度档的设定温度；落在所有区间之外时返回 ``nan``。
        """
        for t_set, step_begin, step_end in lladder:
            if step_begin <= step <= step_end:
                return t_set
        return float('nan')

    def collect_md_frames(self, iround: int = 0, lpotential: list = None,
                          stride: int = 1, n_max_frame: int = None,
                          dict_md_params: dict = None, tag_prefix: str = 'MD',
                          if_skip_if_done: bool = True) -> Path:
        """把本轮所有 MD 轨迹的帧汇成**一份** n2p2 数据集，供 committee 逐帧对比。

        每帧的 ``energy`` / 力写的是**产生它的那条势函数**的 LAMMPS 预测值，只作占位与
        溯源：committee 比较用的是各势的 ``Ennp_phys``（nnp-dataset 现算），真正的
        参考值要等 :meth:`collect_dft` 的 VASP 静态结果。

        Args:
            iround: MD 扩充轮次。
            lpotential: ``[(run, epoch), ...]``；为 None 时读 potentials 的清单。
            stride: 逐轨迹的抽帧步长（1 = 全取）。
            n_max_frame: 汇总后的帧数上限；超出时在**全局**再做一次均匀抽样。
                帧数直接决定后面每条势的 nnp-dataset 代价，务必按预算设上限。
            dict_md_params: 用于取相列表（需与 :meth:`submit_md` 那轮一致）。
            tag_prefix: 帧的 tag 前缀，最终 tag 形如 ``MD0-fcc``。
            if_skip_if_done: 完成阀门；帧集与索引表都非空则直接返回。

        Returns:
            ``data_md/<iround>/y_frames/input.data.frames``。

        Raises:
            FileNotFoundError: 某个 MD 作业目录缺轨迹或温度阶梯文件。
        """
        os.chdir(self.dir_root)
        lpotential = lpotential if lpotential is not None else self._read_lpotential()
        dict_md_params = {**DICT_MD_PARAMS_DEFAULT, **(dict_md_params or {})}
        lphase = dict_md_params['lphase']

        dir_round = self.dir_data_md / str(iround)
        dir_frames = dir_round / 'y_frames'
        path_frames = dir_frames / 'input.data.frames'
        path_index = dir_frames / 'p_post_frames.txt'
        if self._skip_done_stage([path_frames, path_index],
                                 f'collect_md_frames {dir_frames}',
                                 if_skip_if_done=if_skip_if_done):
            return path_frames

        os.makedirs(dir_frames, exist_ok=True)

        # 先全部读进内存再落盘：需要先知道总帧数才能做全局均匀抽样（n_max_frame）。
        # 单帧只有几十个原子，上万帧的 Atoms 列表在内存里是几百 MB 量级，可接受。
        lrecord = []
        for run, epoch in lpotential:
            dir_job = dir_round / 'y_dir' / self._potential_id(run, epoch)
            for phase in lphase:
                path_dump = dir_job / f'y_md_{phase}.lammpstrj'
                lladder = read_temp_ladder(dir_job / f'y_md_{phase}_temp.txt')
                lframe = read_lammps_dump(path_dump, lele=self.lele)
                for step, atoms in lframe[::max(1, int(stride))]:
                    lrecord.append({'run': run, 'epoch': epoch, 'phase': phase,
                                    'T_set_K': self._temperature_of_step(step, lladder),
                                    'step': step, 'natoms': len(atoms),
                                    'E_md_eV': atoms.get_potential_energy(),
                                    'atoms': atoms,
                                    'src': str(path_dump)})
            print(f"  ✅ {self._potential_id(run, epoch)}: "
                  f"{sum(1 for r in lrecord if r['run'] == run and r['epoch'] == epoch)} frame(s)")

        if not lrecord:
            raise FileNotFoundError(f"❌ No MD frames collected under {dir_round / 'y_dir'}.")

        n_all = len(lrecord)
        if n_max_frame is not None and n_all > int(n_max_frame):
            # 全局均匀抽样而非截断：截断会把最后几条势函数的轨迹整条丢掉，
            # committee 的采样覆盖面就偏了。
            lidx = np.linspace(0, n_all - 1, int(n_max_frame)).round().astype(int)
            lrecord = [lrecord[i] for i in sorted(set(lidx.tolist()))]
            print(f"  ⚠️  {n_all} frame(s) exceed n_max_frame={n_max_frame}; "
                  f"uniformly subsampled to {len(lrecord)}")

        # nnp-dataset 的 structure index 是 input.data 里的 0-based 顺序号
        # （Dataset.cpp:816 countStructures 从 0 起、按文件顺序分发），
        # 所以这里的写出顺序就是后面 energy.comp 的 index，两边不能再各自排序。
        writer = nnpdata()
        path_frames.write_text('', encoding='utf-8')
        n_frame = len(lrecord)
        for i, record in enumerate(lrecord):
            record['index'] = i
            writer.write_from_ase(str(path_frames), record['atoms'],
                                  struct_num=i + 1, full_struct_num=n_frame,
                                  file_name=record['src'],
                                  tag=f"{tag_prefix}{iround}-{record['phase']}",
                                  append=True)

        df = pd.DataFrame([{k: r[k] for k in _COLS_FRAMES} for r in lrecord])[_COLS_FRAMES]
        self._write_table(path_index,
                          ['# MD frames collected for committee comparison',
                           f'# round      {iround}',
                           f'# stride     {stride}   n_max_frame {n_max_frame}',
                           f'# phases     {lphase}',
                           '# index      0-based structure index in input.data.frames '
                           '(= nnp-dataset energy.comp index)',
                           '# E_md_eV    total potential energy from the *generating* potential '
                           '(placeholder only)'],
                          df, float_format='16.8f')
        print(f"📊 collect_md_frames: {n_frame} frame(s) from {len(lpotential)} potential(s) "
              f"-> {path_frames}")
        return path_frames


    # ===== 阶段 iii-b：每条势对整份帧集跑 nnp-dataset =====
    def submit_compare(self, iround: int = 0, lpotential: list = None,
                       if_sbatch: bool = False, dict_args_to_submit: dict = None,
                       if_force_rebuild: bool = False, if_skip_if_done: bool = True) -> Path:
        """准备并可选提交 committee 一致性对比作业（每条势函数一个 nnp-dataset）。

        每个作业目录是一份自包含的 nnp 运行环境：势函数三件套 + **整份**
        ``input.data``（= 本轮全部 MD 帧）。``nnp-dataset 0`` 逐结构写
        ``energy.comp.%04d``（index / N / Eref_phys / Ennp_phys），一次跑完所有帧。

        Note:
            ``0`` 是 shuffle 开关，必须为 0：开了洗牌后各 rank 拿到的结构顺序被打乱，
            虽然 index 列仍是全局号，但没有任何好处，只会让并行结果更难核对。
            多 rank 时每个 rank 各写一份 energy.comp，:meth:`collect_compare` 会全部读回
            再按 index 排序。

        Args:
            iround: MD 扩充轮次。
            lpotential: ``[(run, epoch), ...]``；为 None 时读 potentials 的清单。
            if_sbatch: 是否真正提交 Slurm 作业。
            dict_args_to_submit: 传给通用 Slurm 提交器的参数字典。
            if_force_rebuild: 作业目录非空时是否强制清空重建。
            if_skip_if_done: 完成阀门，逐作业生效。

        Returns:
            ``data_md/<iround>/y_compare/``。

        Raises:
            FileNotFoundError: 缺帧集或势函数文件。
        """
        os.chdir(self.dir_root)
        if dict_args_to_submit is None:
            # nnp-dataset 是 MPI 程序；4 核已足够（单帧几十个原子，主要开销是 SF 计算）。
            dict_args_to_submit = {'mode': 'each_subdir', 'chunks': 10,
                                   'module_profile_type': 'zcm6_n2p2_0',
                                   'launcher_type': 'mpirun', 'if_use_my_launcher': False,
                                   'cmd': 'nnp-dataset 0',
                                   'partition': 'amd_512', 'nodes': 1, 'ncores': 4}

        lpotential = lpotential if lpotential is not None else self._read_lpotential()
        dir_round = self.dir_data_md / str(iround)
        path_frames = dir_round / 'y_frames' / 'input.data.frames'
        if not self._is_file_nonempty(path_frames):
            raise FileNotFoundError(f"❌ Missing/empty {path_frames}. Run collect_md_frames first.")

        dir_cmp = dir_round / 'y_compare'
        dir_ydir = dir_cmp / 'y_dir'
        os.makedirs(dir_ydir, exist_ok=True)

        ldir_job = [dir_ydir / self._potential_id(run, epoch) for run, epoch in lpotential]
        ldir_done = ([d for d in ldir_job if any(self._is_file_nonempty(p)
                                                 for p in d.glob('energy.comp*'))]
                     if if_skip_if_done else [])
        if if_skip_if_done and len(ldir_done) == len(ldir_job):
            print(f"⏭️ skip submit_compare {dir_cmp}: all {len(ldir_job)} job dir(s) already carry "
                  f"energy.comp. Pass if_skip_if_done=False to redo them.")
            self.last_jobids = []
            return dir_cmp

        n_kept = 0
        for (run, epoch), dir_job in zip(lpotential, ldir_job):
            if if_skip_if_done and self._is_file_nonempty(dir_job / 'input.data') \
                    and self._is_file_nonempty(dir_job / 'scaling.data'):
                n_kept += 1
                continue
            self._guard_rebuild_dir(dir_job, if_force_rebuild=if_force_rebuild,
                                    label=f'compare job {dir_job.name}')
            os.makedirs(dir_job, exist_ok=True)
            dir_pot_src = self.dir_potentials / run / epoch
            for path_src in [dir_pot_src / 'scaling.data'] \
                    + sorted(dir_pot_src.glob('weights.*.data')):
                shutil.copy(path_src, dir_job / path_src.name)
            self._write_input_nn_energy_only(dir_pot_src / 'input.nn', dir_job / 'input.nn')
            # 帧集逐作业硬拷会把同一份数据复制 N 份（上百 MB × 势数）；nnp-dataset 只读它，
            # 软链接足够，且能一眼看出所有作业吃的是同一份帧集。
            path_link = dir_job / 'input.data'
            if path_link.is_symlink() or path_link.exists():
                path_link.unlink()
            path_link.symlink_to(os.path.relpath(path_frames.resolve(), dir_job.resolve()))
        print(f"Prepared {len(ldir_job) - n_kept} compare job dir(s) under {dir_ydir} "
              f"(kept {n_kept} already-prepared, {len(ldir_done)} of them with results).")

        ldir_todo = [d for d in ldir_job if d not in set(ldir_done)]
        args_to_submit = {**self._clamp_chunks(dict_args_to_submit, len(ldir_todo)),
                          'path_root': str(dir_cmp.resolve()),
                          'dir_root': './y_dir',
                          'if_sbatch': True if if_sbatch else False}
        if ldir_done:
            args_to_submit['lsubdir'] = [d.name for d in ldir_todo]
            print(f"  submit only {len(ldir_todo)} incomplete compare job(s) "
                  f"(skip {len(ldir_done)} already complete)")
        cli_args = self._dict_to_cli_args(args_to_submit)

        self.last_jobids = []
        if if_sbatch:
            before = self._snapshot_jobids()
            os.system(f'pei_slurm_univ_submit.py {cli_args}')
            after = self._snapshot_jobids()
            self.last_jobids = sorted(after - before)
            print(f"  submit_compare observed {len(self.last_jobids)} new Slurm job(s): "
                  f"{self.last_jobids}")
        return dir_cmp


    # ===== 阶段 iii-c / iv：汇总各势预测并按发散度筛结构 =====
    @staticmethod
    def _read_energy_comp(dir_job: Path = None) -> pd.DataFrame:
        """读回单个 nnp-dataset 作业的全部 ``energy.comp.%04d``。

        多 rank 时结构被分到各个 rank、各写一份，必须全部读回再按 index 排序才是
        完整且有序的一份预测；只读 rank 0 会静默丢掉大部分帧。

        Args:
            dir_job: ``y_compare/y_dir/<run>_<epoch>/``。

        Returns:
            列为 ``index`` / ``natoms`` / ``E_nnp_eV`` 的 DataFrame，按 index 升序。

        Raises:
            FileNotFoundError: 目录下没有非空的 energy.comp.*。
        """
        dir_job = Path(dir_job)
        # nnp-dataset 正常跑完会把各 rank 的 energy.comp.%04d 合并成单个 energy.comp 后删掉分片；
        # 作业被砍在合并之前则只剩分片。两种形态都要认，否则「跑完了却读不到」。
        lpath = sorted(p for p in dir_job.glob('energy.comp*') if p.stat().st_size > 0)
        lmerged = [p for p in lpath if p.name == 'energy.comp']
        lpath = lmerged if lmerged else lpath
        if not lpath:
            raise FileNotFoundError(f"❌ No non-empty energy.comp / energy.comp.* in {dir_job}. "
                                    f"nnp-dataset failed? check nnp-dataset.log.*")
        lrow = []
        for path in lpath:
            for line in path.read_text(encoding='utf-8').splitlines():
                if not line.strip() or line.lstrip().startswith('#'):
                    continue
                lfield = line.split()
                # index N Eref_phys Ennp_phys [Eref_int Ennp_int]（nnp-dataset.cpp:277-292）
                lrow.append((int(lfield[0]), int(lfield[1]), float(lfield[3])))
        df = pd.DataFrame(lrow, columns=['index', 'natoms', 'E_nnp_eV'])
        return df.sort_values('index').reset_index(drop=True)

    @staticmethod
    def _write_input_nn_energy_only(path_src: Path = None, path_dst: Path = None) -> None:
        """拷一份 input.nn，同时注释掉 ``use_short_forces``。

        committee 判据只用能量，但 nnp-dataset 见到 ``use_short_forces`` 就会连力一起算
        并落盘 ``forces.comp``：力的反向传播是单点预测里最贵的一步，产物体积也按
        「帧数 × 原子数 × 3」膨胀（万帧量级下每条势就是几十 MB，乘上势数很快上到 10 GB）。
        关掉它对能量预测没有任何影响。

        Args:
            path_src: 势函数库里的 input.nn。
            path_dst: 目标 compare 作业目录里的 input.nn。
        """
        lline = Path(path_src).read_text(encoding='utf-8').splitlines()
        lout = []
        for line in lline:
            if line.split('#', 1)[0].split()[:1] == ['use_short_forces']:
                lout.append('#' + line + '   # disabled by PeiN2p2MD.submit_compare '
                                         '(committee compares energies only)')
            else:
                lout.append(line)
        Path(path_dst).write_text('\n'.join(lout) + '\n', encoding='utf-8')

    def collect_compare(self, iround: int = 0, lpotential: list = None,
                        if_skip_if_done: bool = True) -> Path:
        """逐帧汇总所有势函数的预测能量，写出 committee 一致性表。

        判据就是 chap_4_TL.pdf 的出发点：同一结构如果被训练得好，一组势函数应当给出
        彼此接近的能量；预测越发散，说明该构型在势能面上的表达越不足，越值得补进训练集。

        产物 ``y_compare/y_compare_data_md.txt`` 逐帧一行：帧的溯源信息 + 统计量
        （mean / std / spread）+ **每条势函数一列**的预测能量（eV/atom）。

        Args:
            iround: MD 扩充轮次。
            lpotential: ``[(run, epoch), ...]``；为 None 时读 potentials 的清单。
            if_skip_if_done: 完成阀门。

        Returns:
            ``y_compare/y_compare_data_md.txt``。

        Raises:
            FileNotFoundError: 缺帧索引表，或某个作业没有 energy.comp.*。
            ValueError: 某条势的预测帧数与帧索引表对不上。
        """
        os.chdir(self.dir_root)
        lpotential = lpotential if lpotential is not None else self._read_lpotential()
        dir_round = self.dir_data_md / str(iround)
        dir_cmp = dir_round / 'y_compare'
        path_out = dir_cmp / 'y_compare_data_md.txt'
        if self._skip_done_stage([path_out], f'collect_compare {dir_cmp}',
                                 if_skip_if_done=if_skip_if_done):
            return path_out

        path_index = dir_round / 'y_frames' / 'p_post_frames.txt'
        if not self._is_file_nonempty(path_index):
            raise FileNotFoundError(f"❌ Missing {path_index}. Run collect_md_frames first.")
        df_frame = pd.read_csv(path_index, sep=r'\s+', comment='#')

        df = df_frame[['index', 'phase', 'T_set_K', 'natoms']].copy()
        lcol_pot = []
        for run, epoch in lpotential:
            pot_id = self._potential_id(run, epoch)
            df_pot = self._read_energy_comp(dir_cmp / 'y_dir' / pot_id)
            if len(df_pot) != len(df_frame):
                raise ValueError(f"❌ {pot_id}: nnp-dataset returned {len(df_pot)} structure(s) "
                                 f"but the frame index table has {len(df_frame)}. "
                                 f"input.data mismatch between submit_compare and "
                                 f"collect_md_frames?")
            col = f'E_{pot_id}'
            # 换成 eV/atom：不同帧的原子数可以不同（FCC 32 / HCP 36），总能不可比。
            df[col] = (df_pot['E_nnp_eV'] / df_pot['natoms']).to_numpy()
            lcol_pot.append(col)
            print(f"  ✅ {pot_id}: {len(df_pot)} structure(s)")

        arr = df[lcol_pot].to_numpy(dtype=float)
        df['E_mean_eV/at'] = arr.mean(axis=1)
        # ddof=1：这组势函数是「同一训练配置下的一次随机抽样」，样本标准差才是无偏估计。
        df['E_std_meV/at'] = arr.std(axis=1, ddof=1) * 1e3
        df['E_spread_meV/at'] = (arr.max(axis=1) - arr.min(axis=1)) * 1e3
        df['n_pot'] = len(lcol_pot)
        df = df[_COLS_COMPARE_STAT + lcol_pot]

        self._write_table(path_out,
                          ['# committee consistency of MD frames (n2p2 nnp-dataset, Ennp_phys)',
                           f'# round      {iround}',
                           f'# potentials {len(lcol_pot)} '
                           f'({lcol_pot[0]} .. {lcol_pot[-1]})',
                           '# index      0-based structure index in input.data.frames',
                           '# E_<run>_<epoch>  per-potential prediction, eV/atom',
                           '# E_std / E_spread over the potentials, meV/atom '
                           '(std uses ddof=1)'],
                          df, float_format='16.8f')
        print(f"📊 collect_compare: {len(df)} frame(s) x {len(lcol_pot)} potential(s); "
              f"E_std median {df['E_std_meV/at'].median():.3f} meV/at, "
              f"max {df['E_std_meV/at'].max():.3f} meV/at -> {path_out}")
        return path_out

    def select_md_structures(self, iround: int = 0, n_select_per_group: int = 5,
                             lgroup_by: list = None, std_min: float = None,
                             n_select_total: int = None, tag_prefix: str = 'MD',
                             if_skip_if_done: bool = True) -> Path:
        """按 committee 发散度挑出「训练得不好」的构型，写成待 DFT 标记的清单。

        分组取 top-N 而不是全局取 top-N：发散度与温度强相关（高温帧天然更发散），
        全局排序会让最高的几档温度吃掉全部名额，低温区的欠拟合构型一个都进不来。
        分组口径默认 ``(phase, T_set_K)``，与 chap_4_TL.pdf「每个温度取若干个结构」一致；
        区别是那里按随机取，这里按发散度降序取。

        Args:
            iround: MD 扩充轮次。
            n_select_per_group: 每组取多少个（默认 5，同 LT）。
            lgroup_by: 分组列，默认 ``['phase', 'T_set_K']``；传 ``[]`` 则全局排序。
            std_min: 只保留 ``E_std_meV/at`` 不低于该阈值的帧；None 表示不过滤。
            n_select_total: 全局上限；超出时按发散度降序截断。DFT 标记的代价与它成正比。
            tag_prefix: 写进 input.data.md 的 tag 前缀（tag 形如 ``MD0-fcc``）。
            if_skip_if_done: 完成阀门。

        Returns:
            ``y_compare/input.data.md``。

        Raises:
            FileNotFoundError: 缺 committee 表或帧集。
            ValueError: 过滤后没有任何结构入选。
        """
        os.chdir(self.dir_root)
        lgroup_by = ['phase', 'T_set_K'] if lgroup_by is None else list(lgroup_by)

        dir_round = self.dir_data_md / str(iround)
        dir_cmp = dir_round / 'y_compare'
        path_md = dir_cmp / 'input.data.md'
        path_select = dir_cmp / 'p_post_select.txt'
        if self._skip_done_stage([path_md, path_select], f'select_md_structures {dir_cmp}',
                                 if_skip_if_done=if_skip_if_done):
            return path_md

        path_compare = dir_cmp / 'y_compare_data_md.txt'
        if not self._is_file_nonempty(path_compare):
            raise FileNotFoundError(f"❌ Missing {path_compare}. Run collect_compare first.")
        df = pd.read_csv(path_compare, sep=r'\s+', comment='#')

        n_all = len(df)
        if std_min is not None:
            df = df[df['E_std_meV/at'] >= float(std_min)]
            print(f"  std_min {std_min} meV/at: {len(df)}/{n_all} frame(s) kept")

        if lgroup_by:
            df = (df.sort_values('E_std_meV/at', ascending=False)
                    .groupby(lgroup_by, sort=False, group_keys=False)
                    .head(int(n_select_per_group)))
        df = df.sort_values('E_std_meV/at', ascending=False)
        if n_select_total is not None and len(df) > int(n_select_total):
            print(f"  ⚠️  {len(df)} selected exceed n_select_total={n_select_total}; "
                  f"keeping the most divergent ones")
            df = df.head(int(n_select_total))
        if df.empty:
            raise ValueError(f"❌ No structure selected from {path_compare} "
                             f"(std_min={std_min}, n_select_per_group={n_select_per_group}). "
                             f"Loosen the criteria.")

        # 落盘顺序按 index 升序：input.data.md 与 data_dft 的子目录编号一一对应，
        # 用帧序而不是发散度序，事后回看 DFT 目录时不必再查排名。
        df = df.sort_values('index').reset_index(drop=True)

        # 从帧集里按 index 取结构。帧集是本轮唯一的结构来源，逐条重读比缓存 Atoms 更省心。
        path_frames = dir_round / 'y_frames' / 'input.data.frames'
        if not self._is_file_nonempty(path_frames):
            raise FileNotFoundError(f"❌ Missing {path_frames}. Run collect_md_frames first.")
        data = nnpdata()
        data.load_from_datafile(str(path_frames))

        writer = nnpdata()
        path_md.write_text('', encoding='utf-8')
        lrow = []
        for rank, (_, row) in enumerate(df.iterrows(), start=1):
            idx = int(row['index'])
            tag = f"{tag_prefix}{iround}-{row['phase']}"
            writer.write_from_ase(str(path_md), data.latoms[idx],
                                  struct_num=rank, full_struct_num=len(df),
                                  file_name=f'md_round{iround}_frame{idx}',
                                  tag=tag, append=True)
            lrow.append({'rank': rank, 'index': idx, 'phase': row['phase'],
                         'T_set_K': row['T_set_K'], 'natoms': int(row['natoms']),
                         'E_std_meV/at': row['E_std_meV/at'],
                         'E_spread_meV/at': row['E_spread_meV/at'], 'tag': tag})

        df_out = pd.DataFrame(lrow)[_COLS_SELECT]
        self._write_table(path_select,
                          ['# MD structures selected for DFT labelling (committee divergence)',
                           f'# round             {iround}',
                           f'# group_by          {lgroup_by if lgroup_by else "(global)"}',
                           f'# n_select_per_group {n_select_per_group}   '
                           f'std_min {std_min}   n_select_total {n_select_total}',
                           f'# frames scanned    {n_all}   selected {len(df_out)}',
                           '# rank              1 = most divergent; index = frame index in '
                           'input.data.frames',
                           '# energies in input.data.md are NNP placeholders — they are '
                           'overwritten by DFT in collect_dft'],
                          df_out, float_format='16.8f')
        print(f"📊 select_md_structures: {len(df_out)}/{n_all} frame(s) selected "
              f"(E_std {df_out['E_std_meV/at'].min():.3f} .. "
              f"{df_out['E_std_meV/at'].max():.3f} meV/at) -> {path_md}")
        return path_md


    # ===== 阶段 v-a：DFT 标记 =====
    @classmethod
    def _outcar_is_done(cls, path_outcar: Path = None) -> bool:
        """OUTCAR 存在**且**正常跑完时为 True。

        父类的 :meth:`PeiN2p2._check_outcar_finished` 直接 ``open``，用在「还没跑过的
        目录」上会抛 FileNotFoundError；铺目录阶段恰恰要对尚不存在的 OUTCAR 做判断，
        故这里补一层存在性守卫。

        Args:
            path_outcar: OUTCAR 路径。

        Returns:
            OUTCAR 完整时为 True。
        """
        path_outcar = Path(path_outcar)
        if not path_outcar.is_file():
            return False
        return cls._check_outcar_finished(path_outcar)

    def prepare_dft(self, iround: int = 0, dir_vasp_template: Path = None,
                    lfile_vasp: list = None, if_sbatch: bool = False,
                    dict_args_to_submit: dict = None, if_force_rebuild: bool = False,
                    if_skip_if_done: bool = True) -> Path:
        """为每个入选结构铺一个 VASP 静态计算目录并可选提交。

        MD 帧自己没有 DFT 参考值，``input.data.md`` 里的能量/力只是 NNP 占位值；
        直接把它们喂回训练集等于自蒸馏，势函数只会固化自己的偏差。chap_4_TL.pdf
        Stage I 同样是「DFT 标记之后再加入训练集」。

        Args:
            iround: MD 扩充轮次。
            dir_vasp_template: 提供 INCAR / KPOINTS / POTCAR 的模板目录
                （建议直接指到 ``construct_dataset`` 里某个已跑通的静态算例，保证
                ENCUT / GGA / ISMEAR 等与初始训练集完全同口径）。
            lfile_vasp: 从模板目录逐个拷贝的文件名，默认 ``['INCAR', 'KPOINTS', 'POTCAR']``。
            if_sbatch: 是否真正提交 Slurm 作业。
            dict_args_to_submit: 传给通用 Slurm 提交器的参数字典。
            if_force_rebuild: 作业目录非空时是否强制清空重建。
            if_skip_if_done: 完成阀门，逐目录生效（OUTCAR 已正常结束的不重铺不重投）。

        Returns:
            ``data_dft/<iround>/y_dir/``。

        Raises:
            FileNotFoundError: 缺 input.data.md、模板目录或模板文件。
            ValueError: INCAR 不是静态设置（NSW/IBRION 表明要做离子弛豫）。
        """
        from ase.io import write as ase_write

        os.chdir(self.dir_root)
        lfile_vasp = lfile_vasp or ['INCAR', 'KPOINTS', 'POTCAR']
        if dict_args_to_submit is None:
            dict_args_to_submit = {'preset': 'zcm6_vasp_0', 'mode': 'each_subdir', 'chunks': 10}
        if dir_vasp_template is None:
            raise ValueError("❌ dir_vasp_template must be given (INCAR/KPOINTS/POTCAR source).")
        dir_vasp_template = Path(dir_vasp_template)
        for fn in lfile_vasp:
            if not self._is_file_nonempty(dir_vasp_template / fn):
                raise FileNotFoundError(f"❌ Missing/empty {dir_vasp_template / fn}.")
        self._check_incar_static(dir_vasp_template / 'INCAR')

        path_md = self.dir_data_md / str(iround) / 'y_compare' / 'input.data.md'
        if not self._is_file_nonempty(path_md):
            raise FileNotFoundError(f"❌ Missing {path_md}. Run select_md_structures first.")
        data = nnpdata()
        data.load_from_datafile(str(path_md))

        dir_ydir = self.dir_data_dft / str(iround) / 'y_dir'
        os.makedirs(dir_ydir, exist_ok=True)

        ldir_job = [dir_ydir / f'{i:05d}' for i in range(len(data.latoms))]
        ldir_done = ([d for d in ldir_job if self._outcar_is_done(d / 'OUTCAR')]
                     if if_skip_if_done else [])
        if if_skip_if_done and len(ldir_done) == len(ldir_job):
            print(f"⏭️ skip prepare_dft {dir_ydir}: all {len(ldir_job)} VASP dir(s) already "
                  f"finished. Pass if_skip_if_done=False to redo them.")
            self.last_jobids = []
            return dir_ydir

        n_kept = 0
        for atoms, dir_job in zip(data.latoms, ldir_job):
            if if_skip_if_done and self._is_file_nonempty(dir_job / 'POSCAR'):
                n_kept += 1
                continue
            self._guard_rebuild_dir(dir_job, if_force_rebuild=if_force_rebuild,
                                    label=f'dft job {dir_job.name}')
            os.makedirs(dir_job, exist_ok=True)
            # direct=True + vasp5=True：与 construct_dataset 的其余算例同一份 POSCAR 口径。
            # 不排序原子：input.data 的原子序要与 OUTCAR 的力一一对应，排序会打乱对应关系。
            ase_write(str(dir_job / 'POSCAR'), atoms, format='vasp',
                      direct=True, vasp5=True, sort=False)
            for fn in lfile_vasp:
                shutil.copy(dir_vasp_template / fn, dir_job / fn)
        print(f"Prepared {len(ldir_job) - n_kept} VASP dir(s) under {dir_ydir} "
              f"(kept {n_kept} already-prepared, {len(ldir_done)} of them finished).")

        ldir_todo = [d for d in ldir_job if d not in set(ldir_done)]
        args_to_submit = {**self._clamp_chunks(dict_args_to_submit, len(ldir_todo)),
                          'path_root': str((self.dir_data_dft / str(iround)).resolve()),
                          'dir_root': './y_dir',
                          'if_sbatch': True if if_sbatch else False}
        if ldir_done:
            args_to_submit['lsubdir'] = [d.name for d in ldir_todo]
            print(f"  submit only {len(ldir_todo)} unfinished VASP dir(s) "
                  f"(skip {len(ldir_done)} already finished)")
        cli_args = self._dict_to_cli_args(args_to_submit)

        self.last_jobids = []
        if if_sbatch:
            before = self._snapshot_jobids()
            os.system(f'pei_slurm_univ_submit.py {cli_args}')
            after = self._snapshot_jobids()
            self.last_jobids = sorted(after - before)
            print(f"  prepare_dft observed {len(self.last_jobids)} new Slurm job(s): "
                  f"{self.last_jobids}")
        return dir_ydir

    @staticmethod
    def _check_incar_static(path_incar: Path = None) -> None:
        """确认模板 INCAR 是**静态**设置（NSW=0 且 IBRION=-1）。

        主动学习要的是 MD 帧**原样**的能量与力；一旦模板带着离子弛豫跑，落盘的就是
        弛豫后另一个构型的标签，与 input.data 里的坐标对不上，且这种错配在训练误差里
        几乎看不出来。

        Args:
            path_incar: 模板 INCAR 路径。

        Raises:
            ValueError: NSW / IBRION 表明要做离子弛豫。
        """
        import re as _re
        text = Path(path_incar).read_text(encoding='utf-8')
        dict_tag = {}
        for key in ['NSW', 'IBRION']:
            m = _re.search(rf'(?im)^\s*{key}\s*=\s*(-?\d+)', text)
            dict_tag[key] = int(m.group(1)) if m else None
        if dict_tag['NSW'] not in (0, None) or dict_tag['IBRION'] not in (-1, None):
            raise ValueError(f"❌ {path_incar} is not a static setting "
                             f"(NSW={dict_tag['NSW']}, IBRION={dict_tag['IBRION']}). "
                             f"MD frames must be labelled as-is; ionic relaxation would move "
                             f"the atoms away from the structure stored in input.data.")

    def collect_dft(self, iround: int = 0, sample_type: str = 'train',
                    if_skip_if_done: bool = True) -> Path:
        """读回本轮 VASP 静态结果，汇成 ``input.data.dft``。

        逐目录核对 OUTCAR 是否正常结束；没跑完的目录**跳过并逐个列出**，不静默丢掉——
        缺结构会让本轮的训练集比预期小，而这在学习曲线上看不出来。

        Args:
            iround: MD 扩充轮次。
            sample_type: 写进 ``begin`` 行的 n2p2 集合标记（``'train'`` / ``'test'`` / None）。
            if_skip_if_done: 完成阀门。

        Returns:
            ``data_dft/<iround>/y_dir/input.data.dft``。

        Raises:
            FileNotFoundError: 缺 VASP 目录，或没有任何一个算完。
        """
        os.chdir(self.dir_root)
        dir_ydir = self.dir_data_dft / str(iround) / 'y_dir'
        path_dft = dir_ydir / 'input.data.dft'
        if not dir_ydir.is_dir():
            raise FileNotFoundError(f"❌ Missing {dir_ydir}. Run prepare_dft first.")
        ldir_job = sorted(d for d in dir_ydir.iterdir() if d.is_dir())

        # 完成阀门在这里**不能**只看「文件非空」：上一轮很可能是在部分算例还没跑完时收的，
        # input.data.dft 只装了一部分结构。等剩下的算完再跑本方法时，若按非空跳过，
        # 那批结构就永远进不了训练集，而学习曲线上完全看不出来。故要求条数与目录数相等。
        if if_skip_if_done and self._is_file_nonempty(path_dft) \
                and self._count_structures(path_dft) == len(ldir_job):
            print(f"⏭️ skip collect_dft {dir_ydir}: {len(ldir_job)} structure(s) already collected "
                  f"into {path_dft.name}. Pass if_skip_if_done=False to redo it.")
            return path_dft

        path_select = self.dir_data_md / str(iround) / 'y_compare' / 'p_post_select.txt'
        ltag = None
        if self._is_file_nonempty(path_select):
            df_select = pd.read_csv(path_select, sep=r'\s+', comment='#')
            ltag = df_select['tag'].tolist()

        path_dft.write_text('', encoding='utf-8')
        lfail, n_ok = [], 0
        for dir_job in ldir_job:
            path_outcar = dir_job / 'OUTCAR'
            if not self._outcar_is_done(path_outcar):
                lfail.append(dir_job.name)
                continue
            # 目录名就是 input.data.md 里的 0-based 行号（prepare_dft 用 f'{i:05d}' 命名），
            # 按名字取 tag 而不是按枚举序：少一个目录就会让枚举序整体错位、tag 全部串行。
            i = int(dir_job.name)
            tag = ltag[i] if ltag is not None and i < len(ltag) else f'MD{iround}'
            data = nnpdata()
            data.load_from_outcar(str(path_outcar), index='-1:', tag=tag,
                                  comment_file=f'dft_round{iround}_{dir_job.name}')
            data.write(str(path_dft), append=True, sample_type=sample_type)
            n_ok += 1

        print(f"📊 collect_dft round {iround}: {n_ok}/{len(ldir_job)} OUTCAR(s) usable "
              f"-> {path_dft}")
        if lfail:
            print(f"  ❌ {len(lfail)} unfinished/failed VASP dir(s): {' '.join(lfail)}")
            print(f"  ⚠️  input.data.dft is INCOMPLETE. Fix/rerun those dirs, then call "
                  f"collect_dft({iround}) again (the completion valve will not skip it).")
        if n_ok == 0:
            raise FileNotFoundError(f"❌ No finished OUTCAR under {dir_ydir}.")
        return path_dft


    # ===== 阶段 v-b：回灌训练集 =====
    @staticmethod
    def _count_structures(path_data: Path = None) -> int:
        """数一份 n2p2 input.data 里的结构条数（``begin`` 行数）。

        Args:
            path_data: input.data 路径。

        Returns:
            结构数；文件不存在时为 0。
        """
        path_data = Path(path_data)
        if not path_data.is_file():
            return 0
        return sum(1 for line in path_data.read_text(encoding='utf-8').splitlines()
                   if line.split()[:1] == ['begin'])

    def append_to_train(self, liround: list = None) -> Path:
        """用「不变的初始训练集 + 各轮 DFT 标记结构」重建 ``data/train/input.data``。

        每轮都从 ``input.data.ini`` 重建而不是往 ``input.data`` 上继续追加：追加式写法
        一旦重跑某一轮就会把该轮结构写进去两遍，而重复结构在 n2p2 里既不报错也不易察觉，
        只会悄悄改变各 tag 的有效权重。``input.data.ini`` 在第一次调用时从当时的
        ``input.data`` 快照生成，之后**永不改动**。

        Args:
            liround: 要并入的轮次列表（如 ``[0]``、``[0, 1]``）；为 None 或空则只还原初始集。

        Returns:
            ``data/train/input.data``。

        Raises:
            FileNotFoundError: 初始训练集缺失，或某轮的 input.data.dft 不存在。
        """
        os.chdir(self.dir_root)
        path_train = self.dir_data / 'train' / 'input.data'
        path_ini = self.dir_data / 'train' / 'input.data.ini'

        if not self._is_file_nonempty(path_ini):
            if not self._is_file_nonempty(path_train):
                raise FileNotFoundError(f"❌ Missing both {path_ini} and {path_train}; "
                                        f"copy the initial training set in first.")
            shutil.copy(path_train, path_ini)
            print(f"  📌 snapshot initial training set -> {path_ini} (never modified again)")

        shutil.copy(path_ini, path_train)
        n_ini = self._count_structures(path_train)
        print(f"  restored {n_ini} initial structure(s) from {path_ini}")

        liround = liround or []
        dict_n_round = {}
        for iround in liround:
            path_dft = self.dir_data_dft / str(iround) / 'y_dir' / 'input.data.dft'
            if not self._is_file_nonempty(path_dft):
                raise FileNotFoundError(f"❌ Missing/empty {path_dft}. Run collect_dft({iround}) "
                                        f"first.")
            with open(path_train, 'a', encoding='utf-8') as fd:
                fd.write(path_dft.read_text(encoding='utf-8'))
            dict_n_round[iround] = self._count_structures(path_dft)
            print(f"  ✅ round {iround}: appended {dict_n_round[iround]} DFT-labelled structure(s)")

        n_total = self._count_structures(path_train)
        print(f"📊 append_to_train: {n_ini} initial + {sum(dict_n_round.values())} MD "
              f"= {n_total} structure(s) -> {path_train}")
        return path_train
