SLURM submission API
====================

.. automodule:: mymetal.slurm.submit

``pei_slurm_univ_submit``
-------------------------

.. autofunction:: mymetal.slurm.submit.pei_slurm_univ_submit

Contract
~~~~~~~~

``path_root``
   Absolute ``pathlib.Path``. The function changes the process cwd to this path.

``mode``
   ``parallel``、``each_subdir`` 或 ``single_alloc``。

``dir_root`` / ``lsubdir``
   从 ``path_root / dir_root`` 递归寻找所有 ``y_dir``，再选择其一级 case
   目录；``lsubdir`` 可按 basename 过滤。

``chunks`` / ``chunk_parent_layout``
   正整数并发流数量；layout 为 ``auto``、``shared`` 或 ``per_chunk``。

``module_profile_type`` / ``MODULE_BLOCKS``
   前者必须是后者的 key；dict value 是写入脚本的 module block。

``launcher_type`` / ``cmd``
   launcher 为 ``srun``、``mpirun`` 或 ``none``；``cmd`` 是每个 case 的
   命令。

``partition`` / ``nodes`` / ``ncores``
   Slurm resource fields；nodes/ncores 必须是正整数。它们写入生成脚本，
   dry-run 时不会向 scheduler 验证 partition。

``child_wall_time`` / ``parent_wall_time``
   可选 Slurm duration，例如 ``2-00:00:00``。适用 mode 见
   :doc:`../manual/slurm`。

``if_sbatch``
   默认 ``False``，只生成脚本；``True`` 才调用 sbatch retry wrapper。

返回 ``None``。参数缺失、路径/枚举/resource 非法时 helper 会抛
``SystemExit``。创建脚本失败传播 ``OSError``；``if_sbatch=True`` 时还可能
遇到外部 command/scheduler failure。

最小 dry-run
~~~~~~~~~~~~

.. code-block:: python

   from pathlib import Path
   from mymetal.slurm.submit import pei_slurm_univ_submit

   path_root = Path("/absolute/path/to/slurm-demo")
   pei_slurm_univ_submit(
       path_root=path_root,
       mode="each_subdir",
       dir_root=Path("."),
       chunks=2,
       module_profile_type="none",
       launcher_type="none",
       cmd="echo dry-run",
       partition="debug",
       nodes=1,
       ncores=1,
       if_sbatch=False,
       MODULE_BLOCKS={"none": "# no modules"},
   )

副作用只限于 ``<case>/sub_slurm_univ.sh`` 和 ``path_root/slurm/`` 下的生成
脚本；此例不会调用 ``sbatch``、VASP、LAMMPS 或 n2p2。

Related guides
~~~~~~~~~~~~~~

* :doc:`../manual/slurm`
* :doc:`../reference/scripts`
