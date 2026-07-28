.. _tutorial-outcar-batch:

批量汇总 OUTCAR
===============

:Audience: 需要核对一组 VASP 结果的用户
:Time: 10--15 分钟
:Requires: Python、pjvasp_package、仓库内 tracked OUTCAR
:Runs VASP: No

目标
----

使用 :class:`mymetal.post.newmain.PostTime`、
:class:`mymetal.post.newmain.PostData` 和
:class:`mymetal.post.newmain.PostData2` 读取两个已有 case，得到 convergence、
iteration、elapsed time、energy、stress、pressure、volume 与 maximum force
的确定性 CSV/terminal table。

本教程不调用当前默认路径会失败的 ``post_general()``，也不调用 shell 版
``pei_vasp_univ_post``。

数据与许可边界
--------------

只选择两个 case，避免复制整套约 40 MB OUTCAR：

.. code-block:: text

   mymetal/example/test-post/y_dir/
   ├── 0.997/OUTCAR
   └── 1.000/OUTCAR

原始文件 SHA-256：

.. code-block:: text

   08E01E1C324FB8C2E3CB0FC76B7203FA04FBAA774E7B2E6E484E267C31C705D2  0.997/OUTCAR
   FBE232765B6B22ED5406ACEEAD4CDB803FB094842C402FDB756C2328F3A150FE  1.000/OUTCAR

.. warning::

   仓库当前没有顶层 LICENSE。本教程不把 OUTCAR 复制进 ``docs/``、不提供下载
   archive，也不包含 POTCAR；CI 只读取 clone 中原本已经跟踪的两个文件。维护者
   确认再分发条款后，才能进一步提取并发布独立最小 fixture。

运行完整脚本
------------

.. code-block:: console

   $ python docs/examples/outcar_summary.py \
       mymetal/example/test-post/y_dir \
       -cases 0.997 1.000 \
       -output docs/_build/outcar-summary.csv

.. literalinclude:: ../../examples/outcar_summary.py
   :language: python
   :linenos:

预期输出
--------

.. code-block:: text

   case   status  converged  energy_eV      pressure_kB  volume_A3  fmax_eV_A
   0.997  parsed  True       -23.50620871      5.52000     101.73  0.0000147
   1.000  parsed  True       -23.50776960      0.17000     102.00  0.0000174
   wrote: .../docs/_build/outcar-summary.csv

CSV 还保留 iteration、elapsed time、electronic entropy 和 stress：

.. list-table:: Parsed fixture values
   :header-rows: 1

   * - case
     - converged
     - iteration
     - elapsed (s)
     - energy (eV)
     - EENTRO (eV)
     - stress xx (kB)
     - pressure (kB)
     - volume (Å³)
     - Fmax (eV/Å)
   * - 0.997
     - True
     - 1(11)
     - 800.512
     - -23.50620871
     - 0.00429260
     - 8.25943
     - 5.52
     - 101.73
     - 0.0000146969
   * - 1.000
     - True
     - 1(11)
     - 865.721
     - -23.50776960
     - 0.00364126
     - 0.23318
     - 0.17
     - 102.00
     - 0.0000173781

字段从哪里来
------------

``PostTime``
   查找 ``reached required accuracy``、``Iteration``、``Elapsed time``、
   core group、memory 和最后一个 ``LOOP`` real time。

``PostData``
   读取最后一个 ``energy(sigma->0)``、``EENTRO`` 与六个 ``in kB`` stress
   分量。本页表格只展示 xx；CSV 可按需要继续扩展其余分量。

``PostData2``
   读取 ``volume of cell``、``external pressure``，并从
   ``TOTAL-FORCE (eV/Angst)`` 记录计算 maximum force。

验证方法
--------

* case 数必须等于 ``-cases`` 的数量；
* 两个 tracked fixture 必须得到 ``status=parsed``；
* convergence、energy、pressure、volume 和 Fmax 与上表一致；
* CSV 使用 UTF-8、固定列顺序和确定性 case 顺序；
* 脚本不修改 OUTCAR，也不运行 VASP。

缺少某个 ``OUTCAR`` 时，脚本保留该 case 并输出 ``status=missing``，数值列留空；
不会静默跳过。parser class 的 ``read_OUTCAR(job_list=...)`` 当前必须显式传入
``job_list``；默认 ``None`` 不能直接迭代。

常见错误
--------

``TypeError: 'NoneType' object is not iterable``
   不要直接调用 ``read_OUTCAR()``；传入明确的 case list，或复用本教程脚本。

``status=missing``
   检查目录是否严格为 ``<path_data>/<case>/OUTCAR``，以及大小写是否一致。

字段为空或数值异常
   当前 parser 依赖 VASP 文本标记。先在原文件中核对相应 marker，再决定是否是
   VASP 版本差异或未完成输出，不要为凑表格伪造值。

下一步
------

在 :doc:`../manual/vasp` 查看完整目录 lifecycle；需要批量脚本时查
:doc:`../reference/scripts`。HPC 上提交前先阅读 :doc:`../manual/slurm`。

Related API
-----------

* :class:`mymetal.post.newmain.PostTime`
* :class:`mymetal.post.newmain.PostData`
* :class:`mymetal.post.newmain.PostData2`
