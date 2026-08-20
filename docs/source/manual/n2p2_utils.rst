n2p2_utils 神经网络势函数工具
==============================

``n2p2_utils`` 是 n2p2 神经网络势函数训练的辅助脚本集合，覆盖从 VASP 数据
准备、对称函数生成、active learning 特征选择到训练和后处理的完整流程。

.. note::

   部分脚本依赖外部库 (``StructureForge``, ``skcosmo``, ``yaya_quant``)，
   可能需要额外安装。``n2p2_universal/`` 下的 runner 脚本依赖集群上的
   n2p2 可执行文件 (``nnp-train``, ``nnp-scaling``)。

.. contents:: 本页内容
   :local:
   :depth: 2

完整工作流
----------

.. code-block:: text

   VASP OUTCAR
       │
       ▼
   data_read.py / data_input.py    ← 数据准备
       │
       ▼
   input.data + energy.npy + forces_cfg
       │
       ▼
   sfs_gen_basic_SF.py              ← 对称函数生成
       │
       ▼
   sflist → input.nn
       │
       ▼
   active_sf_0 (sub_cal → collect → select)   ← active learning SF 选择
       │
       ▼
   input.nn (精选 SF)
       │
       ▼
   train_train.bash (nnp-train)     ← 训练
       │
       ▼
   learning-curve.out + weights.*.out
       │
       ▼
   train_get_result.bash + train_final_result.bash   ← 结果提取

数据准备
--------

data_input.py
~~~~~~~~~~~~~

从 VASP 计算目录收集结构，写入 n2p2 ``input.data`` 格式。使用
``StructureForge.InputTools.VASPReader`` 读取每个目录的 ``CONTCAR`` 和
``OUTCAR``。

.. code-block:: python

   # 修改 data 变量指向 VASP 计算目录
   data = "/path/to/vasp/calculations/"
   # 输出 input.data

data_read.py
~~~~~~~~~~~~

读取 VASP 结构和参考力/能量，保存为 ``energy.npy`` 和 ``forces_cfg``。
使用 ``myvasp.vasp_func`` 读取结构，``StructureForge`` 读取力和能量。

对称函数生成
------------

sfs_gen_basic_SF.py
~~~~~~~~~~~~~~~~~~~

生成 n2p2 对称函数 (symmetry function) 参数。支持 G2 (径向) 和 G4 (角向) 类型。

.. code-block:: bash

   python sfs_gen_basic_SF.py -G G2 -e Mg -c 3.3 -n 5
   python sfs_gen_basic_SF.py -G G4 -e Mg -c 4.6 -n 5 -z 1,4,16

.. list-table::
   :header-rows: 1
   :widths: 15 85

   * - 参数
     - 说明
   * - ``-G``
     - 对称函数类型：``G2`` (径向) 或 ``G4`` (角向)
   * - ``-e``
     - 元素符号 (如 ``Mg``, ``Au``)
   * - ``-c``
     - 截断半径 (Å)
   * - ``-n``
     - 点数 (生成 ``2N+1`` 个对称函数)
   * - ``-z``
     - zeta 列表 (默认 ``[1, 4, 16]``，仅 G4)

G2 和 G4 的区别：

- **G2** (径向)：只考虑两个原子的距离，参数为 ``eta`` 和 ``r_shift``
- **G4** (角向)：考虑三个原子的夹角，额外参数为 ``lambda``、``zeta``、``eta``

sfs_gen_gen_sf.bash
~~~~~~~~~~~~~~~~~~~~

批量调用 ``sfs_gen_basic_SF.py`` 生成不同元素、截断半径的对称函数列表，
输出到 ``sflist`` 文件。

.. code-block:: bash

   $SFGEN -G G2 -e Mg -c 3.3 -n 5    > sflist
   $SFGEN -G G2 -e Mg -c 4.6 -n 5  >> sflist
   $SFGEN -G G4 -e Mg -c 4.6 -n 5  >> sflist

Active Learning SF 选择
-----------------------

当对称函数候选过多时，用 active learning 选择最重要的子集。

active_sf_0_sub_cal.py
~~~~~~~~~~~~~~~~~~~~~~

分批计算对称函数子集。从 ``input.nn.all`` 分离所有 ``symfunction_short``
行到 ``SFs_all.dat``，生成不含 SF 的 ``input.nn.no``，然后逐批写入并计算。

active_sf_0_collect_sf.py
~~~~~~~~~~~~~~~~~~~~~~~~~

收集所有对称函数的计算结果。从 ``work_dir/*/function.data`` 读取每个 SF 的
值，提取 ``feat_atom`` (每原子特征) 和 ``feat_av`` (每帧平均特征)。

active_sf_0_select_feat.py
~~~~~~~~~~~~~~~~~~~~~~~~~~

用 CUR 或 FPS 算法 (``skcosmo``) 选择最重要的 N 个对称函数。

.. code-block:: bash

   python active_sf_0_select_feat.py 48    # 选择 48 个最重要的 SF

选择策略：

- **CUR_feat_av**：基于每帧平均特征做 CUR 分解
- **CUR_feat_atom**：按元素分组，每元素选 N/len(elements) 个 SF
- 过滤掉 0 值比例超过 ``max_0`` (5%) 的 SF

active_sf_0_select.bash
~~~~~~~~~~~~~~~~~~~~~~~

active learning SF 选择的入口脚本：

.. code-block:: bash

   $mypython collect_sf.py       # 收集
   $mypython select_feat.py 48   # 选择 48 个
   cp input.nn ../train/          # 复制到训练目录

训练
----

train_train.bash
~~~~~~~~~~~~~~~~

n2p2 训练入口。使用 MPI 并行运行 ``nnp-train``。

.. code-block:: bash

   mpirun -np 16 nnp-train

train_get_result.bash
~~~~~~~~~~~~~~~~~~~~~

获取训练结果：复制 ``input.nn``、``scaling.data``、``learning-curve.out``
到 ``result/post/``，准备后处理。

train_final_result.bash
~~~~~~~~~~~~~~~~~~~~~~~

从训练结果中提取最终势函数：读取最佳 epoch，提取 ``weights.*.out`` 和
``scaling.data`` 到 ``nnp-data`` 目录。

n2p2_universal
--------------

pei_n2p2_univ_run
~~~~~~~~~~~~~~~~~

n2p2 pipeline runner，在 Slurm 分配内运行 ``nnp-scaling → nnp-train``。

.. code-block:: bash

   pei_n2p2_univ_run <job_dir> [scaling_bins] [done_file] [done_marker]

退出码约定：

.. list-table::
   :header-rows: 1
   :widths: 15 85

   * - 退出码
     - 含义
   * - ``0``
     - pipeline 完成
   * - ``10``
     - 已完成，跳过 (检测到哨兵文件)
   * - ``1``
     - 失败或输入缺失

n2p2 没有原生完成标记，runner 写入自己的哨兵文件 (``done_file``) 用于跳过/续算。
``nnp-scaling`` 在 ``scaling.data`` 已存在时跳过；``nnp-train`` 不支持断点续传。

pei_n2p2_univ_load_env
~~~~~~~~~~~~~~~~~~~~~~

n2p2 运行环境加载脚本 (**source 而非执行**)。加载 eigen/gsl/openmpi module，
将 ``nnp-*`` 加入 ``PATH``。

.. code-block:: bash

   source pei_n2p2_univ_load_env

可通过环境变量覆盖默认路径：

- ``N2P2_MODULE_FILE``：module 文件路径
- ``N2P2_MODULES``：要加载的 module 列表
- ``N2P2_BIN_DIR``：n2p2 可执行文件目录

pei_n2p2_univ_clean_train
~~~~~~~~~~~~~~~~~~~~~~~~~

训练目录清理脚本。两阶段清理：

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - 阶段
     - 操作
   * - 阶段 1 (``y_dir/<job>/``)
     - 裁剪 epoch 快照：``neuron-stats``/``trainpoints``/``testpoints``/
       ``trainforces``/``testforces`` 每 50 步留 1 个。
       ``weights.<elem>.<epoch>.out`` 默认不动 (势函数本体)。
   * - 阶段 2 (``y_post/*/properties/y_epoch_scan/y_dir/``)
     - 回收逐 epoch 的 LAMMPS 计算目录 (势函数和汇总数据已提取)。
       删之前通过 "epoch 覆盖校验"。

选项：

- ``--skip-train``：跳过阶段 1
- ``--skip-post``：跳过阶段 2
- ``--include-weights``：对 weights 也按元素分组各留最终 epoch

.. warning::

   此脚本会删除文件。运行前先在非生产目录验证，确认 ``--dry-run`` 输出
   符合预期。

示例脚本
--------

.. literalinclude:: ../../examples/n2p2_utils_overview.py
   :language: python
   :caption: n2p2_utils 子包功能总览示例

.. seealso::

   - :doc:`n2p2` — n2p2 简版说明
   - :doc:`slurm_utils` — Slurm 提交入口 (n2p2 preset)
   - :doc:`lmp_utils` — LAMMPS 势函数计算
   - :doc:`workflows` — Workflow 总览
   - :doc:`../api/ml` — mymetal.ml.n2p2 Python API
