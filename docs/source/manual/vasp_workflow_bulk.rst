vasp_workflow_bulk bulk 目录型 workflow
==========================================

``vasp_utils/vasp_workflow_bulk`` 是 VASP bulk 计算的目录型 workflow 脚本集。
所有 workflow 统一从 ``y_full_relax`` 出发，准备各类后处理目录（EOS、stretch、
NEB、convergence、Cij、cohesive、DOS/band、surface energy、HOEC、KPAR/NCORE），
只准备目录与输入文件，不 sbatch（``pei_vasp_run_cohesive`` 是例外，见下文）。

.. note::

   这些脚本不注册 ``console_scripts``。请将 ``vasp_utils/vasp_workflow_bulk/``
   加入 ``PATH``，或使用完整路径调用。首次使用应先 ``-h`` 或阅读源码头部注释。
   每个 ``run_*`` 脚本都带 ``confirm_prepare_dir`` 闸：若目标 ``y_*`` 已存在且
   没有终端可交互，直接 ``exit 1``，不会在无人值守时删掉旧结果。

.. contents:: 本页内容
   :local:
   :depth: 2

统一入口：y_full_relax
----------------------

所有 bulk workflow 都从 ``y_full_relax`` 出发：它是一个**已收敛**的完整弛豫结果，
提供 ``CONTCAR``（结构源）、``INCAR``、``KPOINTS``、``POTCAR`` 和 ``sub*`` Slurm
脚本。每个 ``run_*`` 脚本都期望在 ``y_full_relax`` 同级目录运行，把生成的批量
目录写到 ``y_<workflow>/y_dir/<case>``。

.. list-table::
   :header-rows: 1
   :widths: 25 25 50

   * - workflow 目录
     - run 脚本
     - 物理量
   * - ``y_eos``
     - ``pei_vasp_run_eos.py``
     - 方程 of state（V/V0 缩放）
   * - ``y_stretch``
     - ``pei_vasp_run_stretch.py``
     - 平衡晶格常数（应变扫描）
   * - ``y_neb``
     - ``pei_vasp_run_neb``
     - NEB 迁移能垒
   * - ``y_convergence``
     - ``pei_vasp_run_convergence``
     - ENCUT/KPOINTS 收敛
   * - ``y_cij_energy``
     - ``pei_vasp_run_cij_energy``
     - 二阶弹性常数（能量-应变法）
   * - ``y_cohesive``
     - ``pei_vasp_run_cohesive``
     - cohesive energy（scale 扫描）
   * - ``y_dos_band``
     - ``pei_vasp_run_dos_band``
     - DOS + band 结构
   * - ``y_surface_energy``
     - ``pei_vasp_run_surface_energy``
     - surface energy（加真空层）
   * - ``y_hoec_energy``
     - ``pei_vasp_run_hoec_energy``
     - 高阶弹性常数（2/3/4 阶）
   * - ``y_kpar_ncore``
     - ``pei_vasp_run_kpar_ncore``
     - KPAR/NCORE 并行基准

.. figure:: /_static/images/generated/vasp_workflow_bulk_overview.png
   :alt: vasp_workflow_bulk 子包 workflow 总览

   ``vasp_utils/vasp_workflow_bulk`` 子包功能总览：10 个 run 脚本从
   ``y_full_relax`` 出发准备 10 类后处理目录，6 个 plot 脚本由
   ``pei_vasp_plot_all`` 调度。

EOS workflow
------------

``pei_vasp_run_eos.py`` 按**体积比** ``V/V0`` 缩放 cell，生成 EOS（equation of
state）计算目录。每个 cell 沿三个方向按 ``ratio^(1/3)`` 缩放（保持形状），
``ISIF`` 默认 4（弛豫形状）。

参数
~~~~

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - 参数
     - 说明
   * - ``-isif``
     - VASP ISIF tag，2 或 4，默认 4
   * - ``-init``
     - 初始体积比，默认 ``0.94``
   * - ``-final``
     - 终止体积比，默认 ``1.06``
   * - ``-interval``
     - 体积比步长，默认 ``0.02``
   * - ``-ratios``
     - 显式体积比列表（逗号分隔），覆盖 ``-init/-final/-interval``
   * - ``-deleteold``
     - 自动删除已存在目录

.. code-block:: bash

   pei_vasp_run_eos.py                                  # 默认 V/V0 = 0.94 .. 1.06, step 0.02
   pei_vasp_run_eos.py -init 0.90 -final 1.10 -interval 0.01   # 更宽更密
   pei_vasp_run_eos.py -isif 2                          # 固定形状
   pei_vasp_run_eos.py -ratios 0.94,0.96,0.98,1.0,1.02,1.04,1.06

Birch-Murnaghan 拟合
~~~~~~~~~~~~~~~~~~~~~

提交 ``y_eos/y_dir/*`` 完成 VASP 计算后，用 ``E(V)`` 数据做 Birch-Murnaghan
（或 Murnaghan）拟合，得到平衡体积 ``V0``、体弹模量 ``B0`` 和 ``B0'``。注意
``-init/-final/-interval`` 是**体积比**（不是应变），cell 按 ``ratio^(1/3)``
缩放。脚本只准备目录，不 sbatch。

.. note::

   ``pei_vasp_plot_all`` 目前没有 ``-eos`` 选项，EOS 拟合仍需手动完成。

.. seealso::

   :doc:`../tutorials/eos_curve` 用合成 Cu-like 数据演示 Murnaghan 与
   Birch-Murnaghan 拟合流程。

Stretch workflow
----------------

``pei_vasp_run_stretch.py`` 沿指定方向施加应变扫描，找平衡晶格常数。与 EOS
不同，这里 ``-init/-final/-interval`` 是**应变**（``-0.004`` = -0.4%），不是
体积比。

参数
~~~~

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - 参数
     - 说明
   * - ``-type``
     - **必需**。拉伸方向：``xyz``/``xy``/``xz``/``yz``/``x``/``y``/``z``
   * - ``-init``
     - 初始应变，默认 ``-0.004`` （-0.4%）
   * - ``-final``
     - 终止应变，默认 ``0.004`` （+0.4%）
   * - ``-interval``
     - 应变步长，默认 ``0.0005`` （0.05%）
   * - ``-strains``
     - 显式应变列表（逗号分隔），覆盖 ``-init/-final/-interval``
   * - ``-keepvolume``
     - 调整未拉伸方向以保持体积不变
   * - ``-deleteold``
     - 自动删除已存在目录

.. code-block:: bash

   pei_vasp_run_stretch.py -type xyz                    # 各向同性，默认 ±0.4% / 0.05%
   pei_vasp_run_stretch.py -type xy -keepvolume         # 面内双轴，体积守恒
   pei_vasp_run_stretch.py -type z -init -0.02 -final 0.02 -interval 0.005
   pei_vasp_run_stretch.py -type xyz -strains=-0.004,-0.002,0.0,0.002,0.004

.. warning::

   当 ``-strains`` 列表以负号开头时（如 ``-0.004,...``），必须用 ``-strains=...``
   形式（等号紧贴），否则 ``argparse`` 会把负数列表误判为 flag。

``-keepvolume`` 说明
~~~~~~~~~~~~~~~~~~~~

``-keepvolume`` 只对**未拉伸**方向做重缩放以守恒体积。例如 ``-type xy
-keepvolume``：xy 面内施加应变，z 方向自动调整使 ``V`` 不变。这对双轴应变
（biaxial）分析尤其有用。

.. seealso::

   :doc:`../tutorials/biaxial_stretch` 用合成数据演示单轴/双轴应变扫描与
   平衡晶格常数拟合。

NEB workflow
------------

``pei_vasp_run_neb``（bash）从 ``y_full_relax_neb/{ini,fin}/CONTCAR`` 生成
``N`` 个 NEB image 目录。

.. code-block:: bash

   # 先用 dist.pl 量端点距离
   dist.pl y_full_relax_neb/ini/CONTCAR y_full_relax_neb/fin/CONTCAR
   # 3.117

   # 推荐 N_IMAGE = DISTANCE / 0.8
   pei_vasp_run_neb 6

源目录必须是 ``y_full_relax_neb``，且包含 ``ini/CONTCAR`` 与 ``fin/CONTCAR``。
脚本只准备 image 目录，不 sbatch。完成后用 ``pei_vasp_plot_all -neb`` 绘制
迁移能垒。

.. seealso::

   NEB 后处理工具集（``nebbarrier``、``nebef``、``nebmovie``、``neb_plot``）
   见 :doc:`neb_utils`。

Convergence workflow
--------------------

``pei_vasp_run_convergence``（bash）做 ENCUT 与 KPOINTS 收敛测试，源目录为
``y_full_relax``。

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - 参数
     - 说明
   * - ``-encuts <val1> <val2> ...``
     - 不同 ENCUT 值列表
   * - ``-kpoints <kx1> <ky1> <kz1> ...``
     - 不同 K 点网格（必须 3 个一组）

.. code-block:: bash

   pei_vasp_run_convergence -encuts 400 500
   pei_vasp_run_convergence -kpoints 3 3 3 4 4 4
   pei_vasp_run_convergence -encuts 400 500 -kpoints 3 3 3 4 4 4

生成 ``y_convergence_encuts`` 与 ``y_convergence_kpoints`` 两棵子树。完成后用
``pei_vasp_plot_all -convergence`` 绘收敛曲线。

Cij / cohesive / HOEC
---------------------

这三个 workflow 都是弹性常数相关，但方法不同：

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   * - workflow
     - 脚本
     - 方法
   * - Cij（二阶）
     - ``pei_vasp_run_cij_energy``
     - 能量-应变法，7 个应变点（-0.003..0.003）
   * - cohesive
     - ``pei_vasp_run_cohesive``
     - scale 扫描，42 个点（0.60..4.00），自动 sbatch 一个 sequential job
   * - HOEC（2/3/4 阶）
     - ``pei_vasp_run_hoec_energy``
     - Wang-Li 法，cubic 11 modes / hex 20 modes

pei_vasp_run_cij_energy
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   pei_vasp_run_cij_energy                              # 默认应变 -0.003..0.003
   pei_vasp_run_cij_energy -e -0.005 -0.002 0.000 0.002 0.005
   pei_vasp_run_cij_energy -no-submit                   # 只准备目录

默认 7 个应变点：``-0.003 -0.002 -0.001 0.000 0.001 0.002 0.003``。``-no-submit``
只准备目录不提交。

pei_vasp_run_cohesive
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   pei_vasp_run_cohesive                                # 无参数

scale 扫描范围 ``0.60`` 到 ``4.00``：``1.30`` 以下步长 ``0.05``，之后步长
``0.10``，共 42 个静态单点。强制 ``NSW=0``、``IBRION=-1``、``ISIF=2``、
``NELM=200``、``LCHARG=F``。在 ``k=4.00`` 时原子几乎自由，该平台即参考能量。

.. note::

   ``pei_vasp_run_cohesive`` 是本子包中**唯一**会自动 sbatch 的脚本：42 个点
   跑在**一个** sequential job 里，不是 42 个独立作业。

pei_vasp_run_hoec_energy
~~~~~~~~~~~~~~~~~~~~~~~~

``pei_vasp_run_hoec_energy``（python）用 Wang-Li 能量-应变法
（Phys. Rev. B 79, 224102 (2009)）计算 2/3/4 阶弹性常数，cubic 11 modes、
hex 20 modes。它是 ``pei_vasp_run_cij_energy`` 的高阶版本：相同入口、相同
``y_dir/<strain>`` 布局，但变形是张量（需要矩阵开方），所以用 Python 而非 bash。

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - 参数
     - 说明
   * - ``-symmetry``
     - ``cubic`` 或 ``hex``，默认自动检测
   * - ``-emax``
     - 单轴参考模式（cubic A、hex M01）的最大振幅
   * - ``-de``
     - 单轴参考模式的步长
   * - ``-static``
     - ``NSW=0, IBRION=-1``，clamped-ion 常数
   * - ``-small_shear``
     - hex only：原位缩放剪切分量，保持亚稳 HCP branch
   * - ``-no_scale_window``
     - 所有 mode 共用一个窗口（旧行为，会过应变混合 mode）

.. code-block:: bash

   pei_vasp_run_hoec_energy                             # 自动检测对称性
   pei_vasp_run_hoec_energy -symmetry hex
   pei_vasp_run_hoec_energy -emax 0.15 -de 0.0025
   pei_vasp_run_hoec_energy -symmetry hex -static       # HCP clamped-ion
   pei_vasp_run_hoec_energy -symmetry hex -small_shear  # relaxed HCP，亚临界剪切

.. warning::

   四阶常数对 ``KPOINTS``/``ENCUT``/``EDIFF`` 极其敏感（Wang-Li Fig. 1, 2, 6）。
   ``y_full_relax`` 必须是收敛的 conventional cell（密 k 点、高 ENCUT、紧
   ``EDIFF ~ 1e-8``）。脚本复用 ``y_full_relax`` 的 ``KPOINTS``/``ENCUT`` 不变。

DOS / band / surface energy
---------------------------

pei_vasp_run_dos_band
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   pei_vasp_run_dos_band                                # 无参数

从 ``y_full_relax_scf`` 出发（注意是 ``_scf`` 不是 ``y_full_relax``），需要
已收敛的 SCF 输出。前置步骤：

1. ``vaspkit -task 303`` 生成 ``KPATH.in``
2. 把 ``PRIM.VASP`` 复制为 ``POSCAR``，跑 SCF（``ISIF=2, NSW=0``）
3. 运行 ``pei_vasp_run_dos_band``

pei_vasp_run_surface_energy
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   pei_vasp_run_surface_energy                          # 无参数

从 ``y_full_relax`` 出发，在 z 方向加真空层构造 surface model，准备 surface vs
bulk 对比作业。脚本只准备目录，不 sbatch。

KPAR/NCORE 基准
---------------

``pei_vasp_run_kpar_ncore``（python）生成 KPAR/NCORE 并行效率基准：固定 VASP
workload（静态单点、固定 KPOINTS/ENCUT、完全约束 lattice），只变 KPAR 和 NCORE。
``NPAR`` 保持注释状态以免与 NCORE 策略冲突。

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - 参数
     - 说明
   * - ``-pairs``
     - ``KPAR:NCORE`` token 列表，默认 21-pair 矩阵（KPAR ∈ {128,64,32,16,8,4}）
   * - ``-ntasks``
     - 总 rank 数，默认 128；每对须满足 ``KPAR*NCORE | ntasks``
   * - ``-force``
     - 重新生成，不提示

.. code-block:: bash

   pei_vasp_run_kpar_ncore                              # 默认 21-pair 矩阵, 128 ranks
   pei_vasp_run_kpar_ncore -pairs 32:4 16:8            # 只测这两对
   pei_vasp_run_kpar_ncore -ntasks 64 -pairs 16:4 8:8 # 测 64-rank job
   pei_vasp_run_kpar_ncore -pairs 16:4 -force         # 重新生成一对，不提示

.. warning::

   ``-ntasks 64`` 单独使用会 abort：默认矩阵仍是 128-rank 对，而
   ``KPAR*NCORE`` 必须整除 ``ntasks``。降低 ``-ntasks`` 时必须同时收窄
   ``-pairs``。

plot_all 调度器
---------------

``pei_vasp_plot_all``（bash）是所有 plot 脚本的统一入口。在 ``y_full_relax``
同级运行，它用 ``find`` 递归发现匹配的 ``y_*`` 目录，进入每个目录运行对应
plot 脚本（plot 脚本会先调 ``pei_vasp_univ_post`` 刮 OUTCAR）。

.. list-table::
   :header-rows: 1
   :widths: 25 30 45

   * - 选项
     - 目标目录
     - plot 脚本
   * - ``-convergence``
     - ``y_convergence``
     - ``pei_vasp_plot_convergence.py``
   * - ``-stretch``
     - ``y_stretch``
     - ``pei_vasp_plot_stretch.py``
   * - ``-neb``
     - ``y_neb``
     - ``pei_vasp_plot_neb.py``
   * - ``-hoec_energy``
     - ``y_hoec_energy``
     - ``pei_vasp_plot_hoec_energy.py``
   * - ``-kpar_ncore``
     - ``y_kpar_ncore``
     - ``pei_vasp_plot_kpar_ncore``
   * - ``-E_in_1_2_bulk``
     - ``y_E_in_1_2_bulk``
     - 单轴/双轴 bulk 能量

.. code-block:: bash

   pei_vasp_plot_all -convergence      # ENCUT/KPOINTS 收敛曲线
   pei_vasp_plot_all -stretch          # 平衡晶格常数
   pei_vasp_plot_all -hoec_energy      # 解高阶弹性常数

.. note::

   只在作业**完成后**运行 ``pei_vasp_plot_all``：它读 OUTCAR，不提交。

pei_vasp_plot_convergence.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

无参数。从 ``y_convergence`` **内部** 运行（``pei_vasp_plot_all`` 的调用方式），
``chdir ..`` 到 ``y_full_relax`` 同级，读 ``y_convergence_encuts`` /
``y_convergence_kpoints``。

pei_vasp_plot_hoec_energy.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - 参数
     - 说明
   * - ``-fitdeg``
     - 多项式拟合阶数
   * - ``-fitmax``
     - 拟合窗口上界（检查 plateau）
   * - ``-skip_univ_post``
     - 跳过重新刮 OUTCAR（~475 个），复用已有数据
   * - ``-maxorder``
     - 最高拟合阶数

输出到 ``y_hoec_energy/``：

- ``y_post_hoec.txt`` — 常数 + 每 mode 拟合诊断
- ``y_post_hoec.json`` — 机器可读常数
- ``y_post_hoec.pdf`` — 每 mode 能量-应变曲线 + 拟合
- ``y_post_hoec_conv.pdf`` — SOEC/TOEC/FOEC 各分量 vs 拟合窗口（**检查收敛用**）

pei_vasp_plot_kpar_ncore
~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - 参数
     - 说明
   * - ``-skip_univ_post``
     - 复用已刮的 ``y_post_time/data.txt``

输出到 ``y_kpar_ncore/``：

- ``p_post_kpar_ncore.txt`` — 可排序 timing/energy 表 + 最快 KPAR/NCORE
- ``p_post_kpar_ncore.pdf`` — time 与 relative energy vs NCORE，每条 KPAR 一曲线

pei_vasp_plot_neb.py
~~~~~~~~~~~~~~~~~~~~

无参数。从 ``y_neb`` **内部** 运行，``chdir ..`` 到 ``y_full_relax`` 同级，读
image 目录绘制 NEB 迁移能垒。

pei_vasp_plot_stretch.py
~~~~~~~~~~~~~~~~~~~~~~~~

无参数。从 ``y_stretch`` **内部** 运行，``chdir ..`` 到 ``y_full_relax`` 同级，
绘应变-能量曲线并拟合平衡晶格常数。

安全检查清单
------------

1. 先 ``bash -n`` 检查 bash 脚本语法，``python -c "import ast; ast.parse(open('x').read())"`` 检查 python 脚本
2. 在非生产目录 dry-run，检查生成的 INCAR/KPOINTS diff
3. ``pei_vasp_run_cohesive`` 会自动 sbatch——确认 ``sub*sequential.sh`` 正确
4. HOEC 四阶常数对收敛参数极敏感——确认 ``y_full_relax`` 已充分收敛
5. ``pei_vasp_run_kpar_ncore`` 降 ``-ntasks`` 时必须同时收窄 ``-pairs``

.. seealso::

   - :doc:`vasp_universal` — 单目录 runner/post/清理工具
   - :doc:`neb_utils` — NEB 后处理工具集
   - :doc:`workflows` — Workflow 总览
   - :doc:`../tutorials/eos_curve` — EOS 拟合教程
   - :doc:`../tutorials/biaxial_stretch` — 单轴/双轴应变教程
   - :doc:`../reference/scripts` — 完整脚本参考
