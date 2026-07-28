Workflow 总览
=============

组件边界
--------

``mymetal``
   可导入的 Python package：结构、I/O、计算、后处理与 ML 数据工具。

``vasp_utils/``、``slurm_utils/``、``lmp_utils/``、``n2p2_utils/``
   目录型 workflow 与 HPC 脚本；不是 Python package 的 console entry point。

VASP、LAMMPS、n2p2、Slurm、MPI 与 module
   集群提供的外部运行环境，本仓库不安装它们。

``mymetal/example/``
   用于演示和解析验证的 tracked fixture；不是生产计算模板。

推荐生命周期
------------

.. code-block:: text

   structure
      -> source inputs
      -> generated case directories
      -> reviewed submit scripts
      -> queued/running jobs
      -> completed or failed outputs
      -> parsed tables and figures

每一箭头都应有一个可观察检查：原子数、输入 diff、``bash -n``、Slurm state、
``OUTCAR`` convergence 或结果断言。不要把“脚本已生成”写成“作业已提交”，也不要
把“作业结束”自动等同于“电子或离子收敛”。

常见目录
--------

.. code-block:: text

   project/
   |-- y_full_relax/
   |-- y_src/
   |   |-- INCAR
   |   |-- KPOINTS
   |   |-- POTCAR
   |   `-- poscars2/
   `-- y_dir/
       |-- 0.997/
       |-- 1.000/
       `-- 1.003/

``y_src`` 保存源输入，``y_dir/<case>`` 保存独立计算。脚本可能假定这一命名；
运行前以实际源码和 ``-h`` 为准。

按任务选择路径
--------------

* 只构建结构：从 :doc:`../getting_started/au111_slab` 开始。
* 计算已有能量的 surface energy：使用 :doc:`../tutorials/surface_energy`。
* 汇总现有 ``OUTCAR``：使用 :doc:`../tutorials/outcar_batch`。
* 生成或提交 HPC 作业：先读 :doc:`slurm`，再读 :doc:`vasp`。
* 查参数和单位：使用 :doc:`../api`；查脚本位置：使用
  :doc:`../reference/scripts`。
