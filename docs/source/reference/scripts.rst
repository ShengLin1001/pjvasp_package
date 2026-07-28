脚本参考
========

入口规则
--------

仓库没有为这些文件注册 ``console_scripts``。下表中的命令名对应真实文件；
可使用仓库内路径运行，也可将所属目录加入 ``PATH``。首次使用应先查看帮助或
源码头部，并在小型目录上 dry-run。

.. code-block:: bash

   python slurm_utils/slurm_universal/pei_slurm_univ_submit.py -h
   bash vasp_utils/vasp_workflow_bulk/pei_vasp_run_neb -h

Slurm
-----

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - 文件
     - 角色
   * - ``slurm_utils/slurm_universal/pei_slurm_univ_submit.py``
     - 软件无关的 argparse 入口；默认只生成脚本。
   * - ``slurm_utils/slurm_universal/pei_slurm_univ_launch_retry``
     - 对可识别的 Slurm/MPI launcher 失败做有限重试。
   * - ``slurm_utils/slurm_universal/pei_slurm_univ_monitor_error``
     - 扫描运行中作业输出中的致命错误。
   * - ``slurm_utils/slurm_vasp/pei_slurm_univ_vasp_monitor``
     - VASP monitor 的便捷启动器。

VASP
----

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - 目录或文件
     - 角色
   * - ``vasp_utils/vasp_universal/pei_vasp_univ_sbatch``
     - 单个 VASP 目录的 runner；会处理历史输出，提交前必须审查。
   * - ``vasp_utils/vasp_universal/pei_vasp_univ_post``
     - 汇总 ``y_dir`` 计算状态。
   * - ``vasp_utils/vasp_workflow_bulk/``
     - 包含 ``pei_vasp_run_eos.py``、``pei_vasp_run_neb``、
       ``pei_vasp_run_surface_energy`` 等 bulk workflow。
   * - ``vasp_utils/vasp_workflow_planar_defects/``
     - 包含 ``pei_vasp_run_gsfe``、``pei_vasp_run_decohesion`` 和
       ``pei_vasp_run_tilt_bulk_to_slip``。
   * - ``vasp_utils/neb_utils/``
     - 包含 ``pei_vasp_univ_neb_select.py``、
       ``pei_vasp_univ_neb_plot.py`` 和结果提取工具。

LAMMPS 与 n2p2
--------------

``lmp_utils/lmp_universal/pei_lmp_run_properties``
   LAMMPS properties workflow 入口。

``n2p2_utils/n2p2_universal/pei_n2p2_univ_run``
   n2p2 workflow runner。

``n2p2_utils/n2p2_universal/pei_n2p2_univ_load_env``
   n2p2 环境加载脚本。

``n2p2_utils/n2p2_universal/pei_n2p2_univ_clean_train``
   训练目录清理脚本；运行前必须审查目标。

安全边界
--------

* 不带 ``-if_sbatch`` 运行 Slurm CLI，先检查生成的 ``sub_*.sh``。
* 清理、续算和重新提交脚本可能改写输入或输出；不要用 ``-h`` 代替源码审查。
* POTCAR、可执行文件、partition、account、module profile 与 launcher 均由目标
  集群提供，不应硬编码为所有用户通用配置。
