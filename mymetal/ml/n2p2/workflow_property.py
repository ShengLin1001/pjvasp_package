"""workflow_property submodule

属性平衡支线的编排层（交接文档 §5.4 的 ``PeiN2p2Property``）。

重活都在 :mod:`mymetal.post.property`（dump 解析、帧规则、护栏、溯源）；
本模块只做编排：物性作业（父类 ``post_properties``）+ epoch 扫描表
（``post_epoch_scan``）完成后，把**结果齐全**的 epoch 目录逐个导出成
DFT 可标记的候选（POSCAR + 静态标记输入 + 溯源），作为物性检查阶段的
常规产物之一——与逐 epoch 表格/图同一层级，而不是散落的手动步骤。

Functions:
    - PeiN2p2Property: 属性平衡支线调度器（唯一对外类）。
"""

from pathlib import Path
import os

from mymetal.ml.n2p2.workflow import PeiN2p2
from mymetal.post.property import export_property_structures


class PeiN2p2Property(PeiN2p2):
    """n2p2 属性平衡支线调度器。

    继承 :class:`~mymetal.ml.n2p2.workflow.PeiN2p2` 的目录约定、物性作业
    （``post_properties``）、epoch 扫描（``post_epoch_scan``）与 Slurm 工具，
    只新增「属性结构出口」一个环节。训练本身仍调父类方法，不重写。

    Attributes:
        dir_data_prop: 候选结构导出根目录（每个 ``<run>_<epoch>`` 一个子目录）。
    """

    def __init__(self, dir_data_prop: Path = Path('./data_prop'), **kwargs):
        """初始化候选导出目录，其余参数原样透传给父类。

        Args:
            dir_data_prop: 属性候选导出根目录。
            **kwargs: 透传给 :class:`PeiN2p2`（``dir_root`` / ``dir_lmp_utils`` /
                ``path_python`` 等必填项仍由调用方给出）。
        """
        super().__init__(**kwargs)
        self.dir_data_prop = Path(dir_data_prop)
        os.makedirs(self.dir_data_prop, exist_ok=True)

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
            path_out = (Path(dir_out) if dir_out is not None
                        else self.dir_data_prop / f'candidates_{dir_run.name}_{epoch}')
            try:
                export_property_structures(path_root=dir_epoch_scan / epoch,
                                           path_out=path_out, lele=self.lele,
                                           run=dir_run.name, epoch=int(epoch),
                                           round_id=round_id, strain_points=strain_points,
                                           path_vasp_static=path_vasp_static)
                dict_summary['exported'].append((dir_run.name, epoch, str(path_out.resolve())))
            except (Exception, SystemExit) as error:
                print(f"❌ export {dir_run.name}/{epoch}: {error}")
                dict_summary['failed'].append((dir_run.name, epoch, str(error)))
        print(f"================ 📊 post_epoch_export {dir_run.name}: "
              f"{len(dict_summary['exported'])} exported, {len(dict_summary['failed'])} failed, "
              f"{len(dict_summary['skipped'])} skipped")
        return dict_summary
