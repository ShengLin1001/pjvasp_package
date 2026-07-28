精选示例
========

先按任务选择可复现教程；需要研究历史目录约定时，再查看
``mymetal/example/``。这些历史数据是工作示例，不是完整测试套件。

无需外部程序
------------

* :doc:`../getting_started/au111_slab`：生成、写出、读回并绘制 Au(111) slab。
* :doc:`../tutorials/surface_energy`：用仓库已有结构和能量计算表面能。
* :doc:`../tutorials/outcar_batch`：批量解析两个真实 ``OUTCAR``。

HPC workflow
------------

* :doc:`../manual/slurm`：先 dry-run 生成脚本，再决定是否提交。
* :doc:`../manual/vasp`：理解 Python API、workflow 脚本和外部 VASP 的边界。
* :doc:`../manual/lammps`：定位 LAMMPS 模板与 runner。
* :doc:`../manual/n2p2`：定位数据准备与训练脚本。

历史 fixture
------------

``mymetal/example/test-generate-bulk``
   Bulk 和 surface 结构生成。

``mymetal/example/test-surface-energy``
   Bulk/slab 能量与结构目录。

``mymetal/example/test-stretch`` 和 ``mymetal/example/test-post``
   ``y_dir/<case>`` 批量目录和 VASP 输出。

``mymetal/example/test-n2p2-sfparams``
   n2p2 symmetry-function 与 dataset 示例。

教程脚本位于 ``docs/examples/``，其中断言是文档构建的最小 smoke test。
