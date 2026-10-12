"""workflow_property submodule

属性平衡支线的编排层（交接文档 §5.4 的 ``PeiN2p2Property``）。

重活都在 :mod:`mymetal.post.property`（dump 解析、帧规则、护栏、溯源）；
本模块只做编排：物性作业（父类 ``post_properties``）+ epoch 扫描表
（``post_epoch_scan``）完成后，把**结果齐全**的 epoch 目录逐个导出成
DFT 可标记的候选（POSCAR + 静态标记输入 + 溯源），作为物性检查阶段的
常规产物之一——与逐 epoch 表格/图同一层级，而不是散落的手动步骤。

导出之后按轮次走 DFT 单点标记：``select_prop_candidates`` 把各势候选按属性
口径筛选并去重成一份待算清单（``data_prop/<k>/y_select/``），
``prepare_prop_dft`` 再逐条铺 VASP 静态目录（``data_dft/<k>/y_dir/<id>/``）
并按需提交。清单冻结后才发现不该算的条目由 ``exclude_prop_after_select`` 记入
``data_dft/<k>/excluded_after_select.jsonl``，prepare 跳过它们（不改 dft_id 编号）。
算完后 ``collect_prop_dft`` 读回单点（扩胞后写 ``y_dir/input.data.dft`` + 逐条 manifest），
``append_prop_to_train`` 用「只读的起始快照 ``input.data.ini`` + 各轮 input.data.dft」重建
``data/train/input.data``（与 MD 支线 ``collect_dft`` / ``append_to_train`` 同一约定）。

Functions:
    - PeiN2p2Property: 属性平衡支线调度器（唯一对外类）。
"""

from collections import Counter
from pathlib import Path
import json
import os
import shutil
import time

import numpy as np

from mymetal.ml.n2p2.dataset import _choose_inplane_rep, nnpdata
from mymetal.ml.n2p2.workflow import PeiN2p2
# 只借 PeiN2p2MD._count_structures（数 begin 行）：回灌约定与 MD 支线一致，不另写一份。
from mymetal.ml.n2p2.workflow_md import PeiN2p2MD
from mymetal.post.property import (DICT_OUTCAR_STATUS_NOTE, check_static_template,
                                   export_property_structures, find_fault_lost_selected,
                                   get_file_hash, get_static_outcar_status,
                                   read_static_dft_frame, select_property_candidates)

# 扩胞前后每原子能量 / 最大力必须不变（E×M、力平铺）；只容浮点舍入。
DTOL_RESIZE_EV = 1e-9


class PeiN2p2Property(PeiN2p2):
    """n2p2 属性平衡支线调度器。

    继承 :class:`~mymetal.ml.n2p2.workflow.PeiN2p2` 的目录约定、物性作业
    （``post_properties``）、epoch 扫描（``post_epoch_scan``）与 Slurm 工具，
    只新增「属性结构出口」一个环节。训练本身仍调父类方法，不重写。

    Attributes:
        dir_data_prop: 属性候选根目录（一级子目录 = 轮次 k，其下每条势一个
            ``candidates_<run>_<epoch>/`` 与去重清单 ``y_select/``）。
        dir_data_dft: 候选 DFT 单点根目录（一级子目录 = 轮次 k）。
    """

    def __init__(self, dir_data_prop: Path = Path('./data_prop'),
                 dir_data_dft: Path = Path('./data_dft'), **kwargs):
        """初始化候选导出与 DFT 标记目录，其余参数原样透传给父类。

        Args:
            dir_data_prop: 属性候选根目录。
            dir_data_dft: 候选 DFT 单点根目录。
            **kwargs: 透传给 :class:`PeiN2p2`（``dir_root`` / ``dir_lmp_utils`` /
                ``path_python`` 等必填项仍由调用方给出）。
        """
        super().__init__(**kwargs)
        self.dir_data_prop = Path(dir_data_prop)
        self.dir_data_dft = Path(dir_data_dft)
        os.makedirs(self.dir_data_prop, exist_ok=True)
        os.makedirs(self.dir_data_dft, exist_ok=True)

    @classmethod
    def _select_export_epochs(cls, dir_epoch_scan: Path = None, lepoch: list = None,
                              if_planar: bool = True) -> tuple:
        """挑出物性**结果齐全**、可导出的 epoch 目录。

        与 :meth:`PeiN2p2._epoch_dir_is_done` 同一个判据：``if_planar`` 决定
        检查三类还是六类结果文件，必须与物性作业的 ``-if_planar`` 一致。

        Args:
            dir_epoch_scan: ``properties/y_epoch_scan/y_dir/``。
            lepoch: 只导出点名 epoch；None 时扫全部数字目录。
            if_planar: 完成判定是否要求平面缺陷三段结果。

        Returns:
            tuple: ``(ldone, lskip)``——可导出与被跳过的 epoch 目录名列表（升序）。
        """
        dir_epoch_scan = Path(dir_epoch_scan)
        if not dir_epoch_scan.is_dir():
            raise FileNotFoundError(f"❌ Missing epoch scan root: {dir_epoch_scan}")
        lepoch_avail = sorted((d.name for d in dir_epoch_scan.iterdir()
                               if d.is_dir() and d.name.isdigit()), key=int)
        if lepoch is not None:
            lwant = sorted({str(e).zfill(6) for e in lepoch}, key=int)
            lmiss = [e for e in lwant if e not in lepoch_avail]
            if lmiss:
                raise FileNotFoundError(f"❌ Epoch dir(s) not found under {dir_epoch_scan}: {lmiss}")
            lepoch_avail = lwant
        ldone = [e for e in lepoch_avail
                 if cls._epoch_dir_is_done(dir_epoch_scan / e, if_planar=if_planar)]
        lskip = [e for e in lepoch_avail if e not in ldone]
        return ldone, lskip

    def post_epoch_export(self, dir_run: Path = Path('./train/y_n2p2_train/y_dir/001'),
                          lepoch: list = None, if_planar: bool = True,
                          dir_out: Path = None, path_vasp_static: Path = None,
                          round_id: str = None, strain_points: int = 5) -> dict:
        """把已完成物性的 epoch 目录导出成 DFT 候选结构。

        逐 epoch 调 :func:`mymetal.post.property.export_property_structures`；
        单个 epoch 导出失败（如某族属性结果损坏）只记录不拖垮整批。
        只读物性产物，物性目录不做任何写操作。

        Args:
            dir_run: n2p2 训练运行目录（与 ``post_properties`` 同一口径）。
            lepoch: 只导出点名 epoch；None 时扫该 run 全部结果齐全的 epoch。
            if_planar: 完成判定是否要求平面缺陷三段结果，须与物性作业一致。
            dir_out: 导出根目录；None 时用 ``self.dir_data_prop``，每 epoch 写入
                其下 ``candidates_<run>_<epoch>/``。
            path_vasp_static: 静态标记模板（INCAR/KPOINTS/POTCAR）目录；None 时
                只导出 POSCAR 与溯源、不复制 VASP 输入（绝不提交计算）。
            round_id: 写入溯源的训练轮次标识。
            strain_points: 每个弹性模式选的对称应变点数（奇数；0 = 全扫描）。

        Returns:
            dict: ``{'exported': [(run, epoch, dir_out), ...], 'failed': [...], 'skipped': [...]}``。

        Raises:
            FileNotFoundError: epoch 扫描根目录或点名的 epoch 不存在。
        """
        os.chdir(self.dir_root)
        dir_run = Path(dir_run)
        dir_epoch_scan = (dir_run.parent.parent / 'y_post' / dir_run.name
                          / 'properties' / 'y_epoch_scan' / 'y_dir')
        ldone, lskip = self._select_export_epochs(dir_epoch_scan, lepoch, if_planar)
        for epoch in lskip:
            print(f"⏭️ skip epoch {epoch}: planar-complete results missing "
                  f"(run the property suite with -if_planar first)")
        dict_summary = {'exported': [], 'failed': [], 'skipped': lskip}
        for epoch in ldone:
            # dir_out 是导出**根**：多个 epoch 共用同一个输出目录会撞 identity 守卫。
            path_out = ((Path(dir_out) if dir_out is not None else self.dir_data_prop)
                        / f'candidates_{dir_run.name}_{epoch}').resolve()
            try:
                # 0 = 全扫描，与 CLI pei_lmp_export_properties 同一约定；库函数只认 None。
                export_property_structures(path_root=(dir_epoch_scan / epoch).resolve(),
                                           path_out=path_out, lele=self.lele,
                                           run=dir_run.name, epoch=int(epoch),
                                           round_id=round_id,
                                           strain_points=strain_points or None,
                                           path_vasp_static=path_vasp_static)
                dict_summary['exported'].append((dir_run.name, epoch, str(path_out.resolve())))
            except (Exception, SystemExit) as error:
                print(f"❌ export {dir_run.name}/{epoch}: {error}")
                dict_summary['failed'].append((dir_run.name, epoch, str(error)))
        print(f"================ 📊 post_epoch_export {dir_run.name}: "
              f"{len(dict_summary['exported'])} exported, {len(dict_summary['failed'])} failed, "
              f"{len(dict_summary['skipped'])} skipped")
        return dict_summary

    # ===== 候选 DFT 单点标记：选帧去重 -> 铺 VASP 静态目录 -> 提交 =====
    def select_prop_candidates(self, iround: int = 0, lexport: list = None,
                               ldefer_model: list = None) -> Path:
        """把本轮各势导出的候选筛选、去重成一份 DFT 待算清单。

        规则全在 :func:`mymetal.post.property.select_property_candidates`：暂缓族只扣
        弛豫帧、伪深谷弛豫帧改由同族未弛豫帧代替、回到完整晶格的弛豫帧与「错排矢量
        即晶格平移」的未弛豫帧不送（与 reference 等价），其余未弛豫帧一律送；
        几何相同的帧只算一次但保留全部来源。本方法只定目录口径。

        Args:
            iround: 属性轮次 k（模型 Mk 产生的候选）。
            lexport: 导出目录列表；None 时取 ``data_prop/<k>/candidates_*``（按名排序）。
            ldefer_model: 弛豫结果暂缓的族名；None 时用库默认 ``LDEFER_MODEL``。

        Returns:
            ``data_prop/<k>/y_select/``。

        Raises:
            FileNotFoundError: 本轮没有任何导出目录。
        """
        os.chdir(self.dir_root)
        dir_round = self.dir_data_prop / str(iround)
        if lexport is None:
            lexport = sorted(d for d in dir_round.glob('candidates_*') if d.is_dir())
        if not lexport:
            raise FileNotFoundError(f"❌ No candidates_* export under {dir_round}.")
        dir_select = dir_round / 'y_select'
        print(f"================ ▶️ select_prop_candidates: round {iround}, {len(lexport)} export(s)")
        for path_export in lexport:
            print(f"  📁 {path_export}")
        select_property_candidates(lpath_export=[Path(d).resolve() for d in lexport],
                                   path_out=dir_select.resolve(), ldefer_model=ldefer_model,
                                   tag_prefix=f'P{iround}')
        return dir_select

    def exclude_prop_after_select(self, iround: int = 0) -> Path:
        """把冻结清单里「层错面消失」的弛豫帧记入排除表，prepare 跳过它们。

        round 0 的清单选定在 fault_lost 护栏之前，VASP 目录已按 dft_id 铺好，重选会改编号；
        故不动 ``select.jsonl``，只重写 ``data_dft/<k>/excluded_after_select.jsonl``（由
        :func:`mymetal.post.property.find_fault_lost_selected` 从清单与导出确定性生成）。
        已铺目录及其输出原样保留、不删；同族未弛豫帧本就在清单里，不用补。

        Args:
            iround: 属性轮次 k。

        Returns:
            排除表 ``data_dft/<k>/excluded_after_select.jsonl``。

        Raises:
            FileNotFoundError: 本轮还没有 prepare 过（缺 ``data_dft/<k>/select.jsonl``）。
        """
        os.chdir(self.dir_root)
        dir_dft_round = self.dir_data_dft / str(iround)
        path_select = dir_dft_round / 'select.jsonl'
        if not self._is_file_nonempty(path_select):
            raise FileNotFoundError(f"❌ Missing {path_select}. Run prepare_prop_dft first.")
        print(f"================ ▶️ exclude_prop_after_select: round {iround}")
        ldrop = find_fault_lost_selected(path_select=path_select.resolve())
        path_excluded = dir_dft_round / 'excluded_after_select.jsonl'
        path_excluded.write_text(''.join(json.dumps(item, ensure_ascii=False) + '\n' for item in ldrop),
                                 encoding='utf-8')
        print(f"✅ exclude_prop_after_select round {iround}: {len(ldrop)} item(s) -> {path_excluded}")
        return path_excluded

    def prepare_prop_dft(self, iround: int = 0, dir_vasp_template: Path = None,
                         lfile_vasp: list = None, if_sbatch: bool = False,
                         dict_args_to_submit: dict = None, if_force_rebuild: bool = False,
                         if_skip_if_done: bool = True) -> Path:
        """按待算清单逐条铺 VASP 静态目录，生成作业脚本，按需提交。

        目录名就是清单里的 ``dft_id``；POSCAR 按字节拷自导出候选并核对 SHA256，
        保证标签与登记的几何一一对应。清单本身复制一份到 ``data_dft/<k>/``，
        之后清单若变了（id 会重排）直接拒绝，避免把旧目录的结果对到新几何上。
        ``if_sbatch=False`` 时仍调用提交器做 dry run，生成的作业脚本可先审阅。

        Args:
            iround: 属性轮次 k。
            dir_vasp_template: 提供 INCAR / KPOINTS / POTCAR 的静态模板目录
                （与训练集同口径，NSW=0、IBRION=-1）。
            lfile_vasp: 从模板逐个拷贝的文件名，默认 ``['INCAR', 'KPOINTS', 'POTCAR']``。
            if_sbatch: 是否真正提交；False 只生成脚本。
            dict_args_to_submit: 传给 ``pei_slurm_univ_submit.py`` 的参数字典。
            if_force_rebuild: 作业目录非空且需重铺时是否清空重建。
            if_skip_if_done: 完成阀门：OUTCAR 正常结束的目录不重铺不重投，
                已铺好的目录（POSCAR 一致）保留。``excluded_after_select.jsonl``
                里的条目始终既不铺也不投。

        Returns:
            ``data_dft/<k>/y_dir/``。

        Raises:
            FileNotFoundError: 缺待算清单或模板文件。
            ValueError: 模板不是静态设置、清单已变、或已有 POSCAR 与清单不符。
            RuntimeError: 提交引擎拒绝（如本轮 ``data_dft/<k>/slurm/`` 还有作业在队，dry run
                也会拒绝，避免改写在跑的 chunk 脚本）或提交失败。
        """
        ### check ###
        os.chdir(self.dir_root)
        lfile_vasp = lfile_vasp or ['INCAR', 'KPOINTS', 'POTCAR']
        if dict_args_to_submit is None:
            dict_args_to_submit = {'preset': 'zcm6_vasp_0', 'mode': 'each_subdir', 'chunks': 10}
        if dir_vasp_template is None:
            raise ValueError("❌ dir_vasp_template must be given (INCAR/KPOINTS/POTCAR source).")
        dir_vasp_template = Path(dir_vasp_template).resolve()
        for fn in lfile_vasp:
            if not self._is_file_nonempty(dir_vasp_template / fn):
                raise FileNotFoundError(f"❌ Missing/empty {dir_vasp_template / fn}.")
        check_static_template(dir_vasp_template)
        path_select = self.dir_data_prop / str(iround) / 'y_select' / 'select.jsonl'
        if not self._is_file_nonempty(path_select):
            raise FileNotFoundError(f"❌ Missing {path_select}. Run select_prop_candidates first.")
        lunique = [json.loads(line) for line in path_select.read_text(encoding='utf-8').splitlines()
                   if line.strip()]
        # to here

        ### prepare ###
        dir_dft_round = self.dir_data_dft / str(iround)
        dir_ydir = dir_dft_round / 'y_dir'
        os.makedirs(dir_ydir, exist_ok=True)
        path_select_copy = dir_dft_round / 'select.jsonl'
        if path_select_copy.is_file():
            if get_file_hash(path_select_copy) != get_file_hash(path_select):
                raise ValueError(f"❌ {path_select} changed since {dir_ydir} was prepared; dft_id "
                                 f"would no longer match the existing dirs. Use a new round dir.")
        else:
            shutil.copy(path_select, path_select_copy)
        path_excluded = dir_dft_round / 'excluded_after_select.jsonl'
        lexcluded_id = ([json.loads(line)['dft_id'] for line in
                         path_excluded.read_text(encoding='utf-8').splitlines() if line.strip()]
                        if path_excluded.is_file() else [])
        # to here

        ### main ###
        ldone, n_new, n_kept = [], 0, 0
        for item in lunique:
            if item['dft_id'] in lexcluded_id:
                continue
            dir_job = dir_ydir / item['dft_id']
            path_poscar = dir_job / 'POSCAR'
            if if_skip_if_done and self._is_file_nonempty(path_poscar):
                if get_file_hash(path_poscar) != item['poscar_sha256']:
                    raise ValueError(f"❌ {path_poscar} differs from the selected geometry "
                                     f"({item['poscar']}).")
                if self._is_file_nonempty(dir_job / 'OUTCAR') \
                        and self._check_outcar_finished(dir_job / 'OUTCAR'):
                    ldone.append(item['dft_id'])
                n_kept += 1
                continue
            self._guard_rebuild_dir(dir_job, if_force_rebuild=if_force_rebuild,
                                    label=f'property dft job {item["dft_id"]}')
            os.makedirs(dir_job, exist_ok=True)
            if get_file_hash(item['poscar']) != item['poscar_sha256']:
                raise ValueError(f"❌ Source POSCAR changed after selection: {item['poscar']}")
            # 字节拷贝而非 ase 重写：保住导出时做过往返校验的那份几何与其 SHA256。
            shutil.copy(item['poscar'], path_poscar)
            for fn in lfile_vasp:
                shutil.copy(dir_vasp_template / fn, dir_job / fn)
            n_new += 1
        print(f"✅ prepare_prop_dft round {iround}: {n_new} new, {n_kept} kept "
              f"({len(ldone)} finished, {len(lexcluded_id)} excluded after select) "
              f"of {len(lunique)} VASP dir(s) under {dir_ydir}")
        for dft_id in lexcluded_id:
            print(f"  ⏭️ {dft_id}: excluded after select ({path_excluded.name}), not submitted")
        ltodo = [item['dft_id'] for item in lunique
                 if item['dft_id'] not in set(ldone) and item['dft_id'] not in lexcluded_id]
        self.last_jobids = []
        if not ltodo:
            print(f"⏭️ skip submit: all {len(lunique)} VASP dir(s) finished or excluded "
                  f"({len(ldone)} finished, {len(lexcluded_id)} excluded).")
            return dir_ydir
        args_to_submit = {**self._clamp_chunks(dict_args_to_submit, len(ltodo)),
                          'path_root': str(dir_dft_round.resolve()),
                          'dir_root': './y_dir',
                          'if_sbatch': True if if_sbatch else False}
        if ldone or lexcluded_id:
            args_to_submit['lsubdir'] = ltodo
            print(f"  submit only {len(ltodo)} unfinished VASP dir(s) "
                  f"(skip {len(ldone)} finished, {len(lexcluded_id)} excluded)")
        cli_args = self._dict_to_cli_args(args_to_submit)
        # dry run 也走提交器：生成的 sub 脚本（节点/核数/启动器）要在真提交前给人审阅。
        if if_sbatch:
            before = self._snapshot_jobids()
            self._run_submit_cli(cli_args, label='prepare_prop_dft')
            after = self._snapshot_jobids()
            self.last_jobids = sorted(after - before)
            print(f"  prepare_prop_dft observed {len(self.last_jobids)} new Slurm job(s): "
                  f"{self.last_jobids}")
        else:
            self._run_submit_cli(cli_args, label='prepare_prop_dft (dry run)')
            print(f"📍 dry run only (if_sbatch=False): review {dir_ydir} and its job scripts, "
                  f"then rerun with if_sbatch=True.")
        return dir_ydir
        # to here

    # ===== DFT 单点读回 -> 扩胞 -> 回灌训练集 =====
    def _read_prop_dft_list(self, iround: int = 0) -> tuple:
        """读本轮冻结清单与排除表，给出应收条目。

        Args:
            iround: 属性轮次 k。

        Returns:
            tuple: ``(litem, lexcluded_id)``——应收条目（清单顺序、已扣排除表）与排除的 dft_id。

        Raises:
            FileNotFoundError: 本轮还没有 prepare 过（缺 ``data_dft/<k>/select.jsonl``）。
        """
        dir_dft_round = self.dir_data_dft / str(iround)
        path_select = dir_dft_round / 'select.jsonl'
        if not self._is_file_nonempty(path_select):
            raise FileNotFoundError(f"❌ Missing {path_select}. Run prepare_prop_dft first.")
        lunique = [json.loads(line) for line in path_select.read_text(encoding='utf-8').splitlines()
                   if line.strip()]
        path_excluded = dir_dft_round / 'excluded_after_select.jsonl'
        lexcluded_id = ([json.loads(line)['dft_id'] for line in
                         path_excluded.read_text(encoding='utf-8').splitlines() if line.strip()]
                        if path_excluded.is_file() else [])
        litem = [item for item in lunique if item['dft_id'] not in lexcluded_id]
        return litem, lexcluded_id

    @staticmethod
    def _get_nnpdata_frame(atoms=None, tag: str = None, comment_file: str = None,
                           path_source: Path = None):
        """把一帧带标签的 Atoms 包成 nnpdata，供 ``adjust_size`` / ``write`` 复用。

        nnpdata 只有从 OUTCAR / input.data 读的入口；属性帧的几何取自 POSCAR（见
        :func:`mymetal.post.property.read_static_dft_frame`），故按 ``load_from_outcar``
        填的同一组字段手工装入，写出格式与 MD 轮次逐字一致。

        Args:
            atoms: 带 SinglePointCalculator 的 Atoms。
            tag: comment 行的 ``tag=``。
            comment_file: comment 行的 ``file=``。
            path_source: 来源 OUTCAR。

        Returns:
            nnpdata: 单结构数据对象。
        """
        data = nnpdata()
        data.latoms = [atoms]
        data.ltags = [tag]
        data.lcomment_files = [comment_file]
        data.lfiles = [str(path_source)]
        data.lstruct_numbers = [1]
        data.lfull_struct_numbers = [1]
        data.lforces = [atoms.get_forces(apply_constraint=False)]
        data.lenergies = [atoms.get_potential_energy()]
        data.lsample_types = [None]
        data.dict = {'latoms': data.latoms, 'ltags': data.ltags,
                     'lcomment_files': data.lcomment_files, 'lfiles': data.lfiles,
                     'lstruct_numbers': data.lstruct_numbers,
                     'lfull_struct_numbers': data.lfull_struct_numbers,
                     'lforces': data.lforces, 'lenergies': data.lenergies,
                     'lsample_types': data.lsample_types}
        return data

    def collect_prop_dft(self, iround: int = 0, sample_type: str = 'train',
                         if_adjust_size: bool = True, size_top: int = 72, size_bottom: int = 36,
                         size_close: int = 48, if_skip_if_done: bool = True) -> Path:
        """读回本轮属性候选的 VASP 单点，扩胞后汇成 ``input.data.dft``。

        按冻结清单 ``data_dft/<k>/select.jsonl`` 逐条读，排除表里的条目逐个打印后跳过。
        判据比 MD 版严：OUTCAR 正常结束、含 ``EDIFF is reached``、恰好一个离子步；不合格的
        逐个列出（不静默丢）。POSCAR 的 SHA256 必须仍等于清单登记值，否则直接报错——
        标签与几何对不上比缺一帧更糟。tag 用清单里的 ``P<k>-<phase>-<group>-<model>``；
        全部来源、扩胞倍数、原/新原子数写进 ``input.data.dft.manifest.jsonl``。

        扩胞沿用 ``nnpdata.adjust_size``（PJ 2026-10-10：先去重再算、合并时扩胞）：面内
        复制、nz=1、E×M、力平铺；逐帧核对每原子能量与最大力不变。

        Args:
            iround: 属性轮次 k。
            sample_type: 写进 ``begin`` 行的 n2p2 集合标记（``'train'`` / ``'test'`` / None）。
            if_adjust_size: 是否扩胞。
            size_top: 扩胞上限原子数（已超出的保持不变）。
            size_bottom: 扩胞下限原子数。
            size_close: 扩胞目标原子数。
            if_skip_if_done: 完成阀门：input.data.dft 与 manifest 条数都等于应收条数才跳过。

        Returns:
            ``data_dft/<k>/y_dir/input.data.dft``。

        Raises:
            FileNotFoundError: 缺清单 / VASP 目录，或没有任何一条可用。
            ValueError: POSCAR 与清单不符、OUTCAR 坐标与 POSCAR 不符、扩胞改变了每原子能量或最大力。
        """
        ### check ###
        os.chdir(self.dir_root)
        dir_ydir = self.dir_data_dft / str(iround) / 'y_dir'
        path_dft = dir_ydir / 'input.data.dft'
        path_manifest = dir_ydir / 'input.data.dft.manifest.jsonl'
        if not dir_ydir.is_dir():
            raise FileNotFoundError(f"❌ Missing {dir_ydir}. Run prepare_prop_dft first.")
        litem, lexcluded_id = self._read_prop_dft_list(iround)
        # 与 MD 版同理：只看「非空」会把部分收齐的旧文件当成完成，补算的帧就永远进不来。
        if if_skip_if_done and PeiN2p2MD._count_structures(path_dft) == len(litem) \
                and self._is_file_nonempty(path_manifest) \
                and len(path_manifest.read_text(encoding='utf-8').splitlines()) == len(litem):
            print(f"⏭️ skip collect_prop_dft {dir_ydir}: {len(litem)} structure(s) already collected "
                  f"into {path_dft.name}. Pass if_skip_if_done=False to redo it.")
            return path_dft
        # to here

        ### main ###
        print(f"================ ▶️ collect_prop_dft: round {iround}, {len(litem)} item(s) "
              f"({len(lexcluded_id)} excluded after select)")
        for dft_id in lexcluded_id:
            print(f"  ⏭️ {dft_id}: excluded after select, not collected")
        ldata, lrecord, lfail = [], [], []
        for item in litem:
            dir_job = dir_ydir / item['dft_id']
            path_outcar = dir_job / 'OUTCAR'
            status = get_static_outcar_status(path_outcar)
            if status != 'ok':
                lfail.append((item['dft_id'], status))
                continue
            if get_file_hash(dir_job / 'POSCAR') != item['poscar_sha256']:
                raise ValueError(f"❌ {dir_job / 'POSCAR'} differs from the selected geometry "
                                 f"({item['poscar']}); its labels no longer match the list.")
            atoms = read_static_dft_frame(dir_job / 'POSCAR', path_outcar)
            n0, energy = len(atoms), atoms.get_potential_energy()
            fmax = np.abs(atoms.get_forces()).max()
            nx, ny = (_choose_inplane_rep(n0, size_top=size_top, size_bottom=size_bottom,
                                          size_close=size_close) if if_adjust_size else (1, 1))
            data = self._get_nnpdata_frame(atoms, tag=item['tag'],
                                           comment_file=f"dft_prop{iround}_{item['dft_id']}",
                                           path_source=path_outcar.resolve())
            if if_adjust_size:
                data.adjust_size(size_top=size_top, size_bottom=size_bottom, size_close=size_close)
            n1 = len(data.latoms[0])
            if n1 != n0 * nx * ny \
                    or abs(data.lenergies[0] / n1 - energy / n0) > DTOL_RESIZE_EV \
                    or abs(np.abs(data.lforces[0]).max() - fmax) > DTOL_RESIZE_EV:
                raise ValueError(f"❌ {item['dft_id']}: resize {n0}->{n1} changed E/atom or max|F|.")
            ldata.append(data)
            lrecord.append(dict(dft_id=item['dft_id'], tag=item['tag'], sample_type=sample_type,
                                natoms_dft=n0, natoms_data=n1, rep=[nx, ny, 1],
                                energy_eV=energy, energy_eV_atom=energy / n0,
                                free_energy_eV=atoms.calc.get_property('free_energy', atoms),
                                fmax_eV_A=float(fmax), poscar_sha256=item['poscar_sha256'],
                                outcar=str(path_outcar.resolve()), sources=item['sources']))
        # 全部读完、核完再落盘：中途报错不会留下半份 input.data.dft。
        path_dft.write_text('', encoding='utf-8')
        for data in ldata:
            data.write(str(path_dft), append=True, sample_type=sample_type)
        path_manifest.write_text(''.join(json.dumps(record, ensure_ascii=False) + '\n'
                                         for record in lrecord), encoding='utf-8')
        # to here

        ### summary ###
        dict_tag = Counter(record['tag'] for record in lrecord)
        dict_resize = Counter((record['natoms_dft'], record['natoms_data']) for record in lrecord)
        print(f"================ 📊 collect_prop_dft round {iround}: {len(lrecord)}/{len(litem)} "
              f"usable -> {path_dft} (+ {path_manifest.name})")
        for tag, count in sorted(dict_tag.items()):
            print(f"  ✅ {tag:<34s} {count:4d}")
        for (n0, n1), count in sorted(dict_resize.items()):
            print(f"  📐 natoms {n0:3d} -> {n1:3d}: {count:4d} structure(s)")
        if lfail:
            print(f"  ❌ {len(lfail)} unusable VASP dir(s):")
            for dft_id, status in lfail:
                print(f"    ❌ {dft_id}: {status} ({DICT_OUTCAR_STATUS_NOTE[status]})")
            print(f"  ⚠️  input.data.dft is INCOMPLETE. Fix/rerun those dirs, then call "
                  f"collect_prop_dft({iround}) again (the completion valve will not skip it).")
        if not lrecord:
            raise FileNotFoundError(f"❌ No usable OUTCAR under {dir_ydir}.")
        return path_dft
        # to here

    @staticmethod
    def _count_sample_types(path_data: Path = None) -> dict:
        """数 input.data 里 ``begin set=train`` / ``begin set=test`` / 裸 ``begin`` 各多少。

        Args:
            path_data: input.data 路径。

        Returns:
            dict: ``{'train': n, 'test': n, 'unset': n}``。
        """
        dict_count = {'train': 0, 'test': 0, 'unset': 0}
        for line in Path(path_data).read_text(encoding='utf-8').splitlines():
            lpart = line.split()
            if lpart[:1] != ['begin']:
                continue
            dict_count[lpart[1].split('=')[1] if len(lpart) > 1 else 'unset'] += 1
        return dict_count

    def append_prop_to_train(self, liround: list = None, path_data_ini_src: Path = None,
                             sha256_ini_src: str = None) -> Path:
        """用「只读起始快照 + 各轮属性 DFT 标签」重建 ``data/train/input.data``。

        与 :meth:`PeiN2p2MD.append_to_train` 同一约定：``input.data.ini`` 一经快照永不改动，
        每次都从它重建而不是在 input.data 上续写（重跑某轮不会把同一批结构写两遍）。
        stage3 的起始集来自上一 stage，``input.data.ini`` 缺失时必须给出来源文件，
        按字节拷贝（保留每个结构的 ``set=train/test``）并把来源与 SHA256 记进
        ``input.data.ini.source.json``；已有快照时来源若给出则必须逐字节一致。
        每轮的 input.data.dft 必须收齐（条数 = manifest 条数 = 清单减排除表）才并入。

        Args:
            liround: 要并入的属性轮次（如 ``[0]``）；None 或空则只还原起始集。
            path_data_ini_src: 起始集来源（如 ``stage2/0/data/train/input.data``）。
            sha256_ini_src: 来源文件应有的 SHA256；给出时先核对，防止上游数据被改过。

        Returns:
            ``data/train/input.data``。

        Raises:
            FileNotFoundError: 缺起始快照且没给来源，或某轮 input.data.dft 缺失。
            ValueError: 来源 SHA 不符、已有快照与来源不一致，或某轮没收齐。
        """
        ### check ###
        os.chdir(self.dir_root)
        dir_train = self.dir_data / 'train'
        path_train = dir_train / 'input.data'
        path_ini = dir_train / 'input.data.ini'
        path_ini_source = dir_train / 'input.data.ini.source.json'
        if path_data_ini_src is not None:
            path_data_ini_src = Path(path_data_ini_src).resolve()
            if not self._is_file_nonempty(path_data_ini_src):
                raise FileNotFoundError(f"❌ Missing/empty {path_data_ini_src}.")
            sha_src = get_file_hash(path_data_ini_src)
            if sha256_ini_src is not None and sha_src != sha256_ini_src:
                raise ValueError(f"❌ {path_data_ini_src} sha256 {sha_src} != expected {sha256_ini_src}.")
            if self._is_file_nonempty(path_ini) and get_file_hash(path_ini) != sha_src:
                raise ValueError(f"❌ {path_ini} differs from {path_data_ini_src}; the start "
                                 f"snapshot must never change.")
        elif not self._is_file_nonempty(path_ini):
            raise FileNotFoundError(f"❌ Missing {path_ini}; pass path_data_ini_src to snapshot "
                                    f"the start dataset.")
        liround = liround or []
        dict_path_dft = {}
        for iround in liround:
            path_dft = self.dir_data_dft / str(iround) / 'y_dir' / 'input.data.dft'
            path_manifest = path_dft.parent / 'input.data.dft.manifest.jsonl'
            if not self._is_file_nonempty(path_dft) or not self._is_file_nonempty(path_manifest):
                raise FileNotFoundError(f"❌ Missing {path_dft} or its manifest. Run "
                                        f"collect_prop_dft({iround}) first.")
            litem, _ = self._read_prop_dft_list(iround)
            n_dft = PeiN2p2MD._count_structures(path_dft)
            n_manifest = len(path_manifest.read_text(encoding='utf-8').splitlines())
            if not n_dft == n_manifest == len(litem):
                raise ValueError(f"❌ round {iround} incomplete: {n_dft} structure(s), {n_manifest} "
                                 f"manifest line(s), {len(litem)} expected. Finish the DFT dirs and "
                                 f"rerun collect_prop_dft({iround}).")
            dict_path_dft[iround] = path_dft
        # to here

        ### prepare ###
        os.makedirs(dir_train, exist_ok=True)
        if not self._is_file_nonempty(path_ini):
            shutil.copy(path_data_ini_src, path_ini)
            dict_source = dict(source=str(path_data_ini_src), sha256=get_file_hash(path_ini),
                               n_structure=PeiN2p2MD._count_structures(path_ini),
                               sample_types=self._count_sample_types(path_ini),
                               created=time.strftime('%Y-%m-%d %H:%M:%S'))
            path_ini_source.write_text(json.dumps(dict_source, ensure_ascii=False, indent=2) + '\n',
                                       encoding='utf-8')
            print(f"  📌 snapshot start dataset {path_data_ini_src} -> {path_ini} "
                  f"(sha256 {dict_source['sha256'][:16]}…, never modified again)")
        # to here

        ### main ###
        shutil.copy(path_ini, path_train)
        n_ini = PeiN2p2MD._count_structures(path_train)
        print(f"  restored {n_ini} start structure(s) from {path_ini} {self._count_sample_types(path_ini)}")
        dict_n_round = {}
        for iround, path_dft in dict_path_dft.items():
            with open(path_train, 'a', encoding='utf-8') as fd:
                fd.write(path_dft.read_text(encoding='utf-8'))
            dict_n_round[iround] = PeiN2p2MD._count_structures(path_dft)
            print(f"  ✅ property round {iround}: appended {dict_n_round[iround]} DFT-labelled "
                  f"structure(s) from {path_dft}")
        n_total = PeiN2p2MD._count_structures(path_train)
        if n_total != n_ini + sum(dict_n_round.values()):
            raise ValueError(f"❌ {path_train} holds {n_total} structure(s), expected "
                             f"{n_ini} + {sum(dict_n_round.values())}.")
        print(f"================ 📊 append_prop_to_train: {n_ini} start + {sum(dict_n_round.values())} "
              f"property = {n_total} structure(s) {self._count_sample_types(path_train)} -> {path_train} "
              f"(sha256 {get_file_hash(path_train)[:16]}…)")
        return path_train
        # to here
