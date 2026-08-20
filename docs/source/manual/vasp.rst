VASP workflow
=============

本页是 VASP workflow 的图文总览：从应用场景、输入文件结构、执行流程到基于集群
真实计算结果的展示。子模块的逐脚本说明见 :doc:`vasp_universal` 与
:doc:`vasp_workflow_bulk`。

.. contents:: 本页内容
   :local:
   :depth: 2

功能分层
--------

``mymetal.io.vasp``
   Python 层读写 POSCAR/CONTCAR；不运行 VASP。

``mymetal.post``
   解析已存在的 ``OUTCAR`` 与批量目录；不判断调度器状态。

``vasp_utils/vasp_universal``
   单目录 runner、monitor、post、续算和清理脚本。

``vasp_utils/vasp_workflow_*``
   EOS、stretch、NEB、surface energy、GSFE 等目录型 workflow。

``myvasp``
   历史 shell 与 VTST-style 工具。可导入的 ``myvasp`` Python helper 来自
   ``myalloy_package`` 。

应用场景
--------

VASP workflow 覆盖以下材料性质计算场景，每个场景对应一个 ``run_*`` 脚本：

.. list-table::
   :header-rows: 1
   :widths: 22 28 50

   * - workflow
     - 脚本
     - 应用场景
   * - EOS
     - ``pei_vasp_run_eos.py``
     - 方程 of state，拟合体弹模量 ``B0`` 和平衡体积 ``V0``
   * - stretch
     - ``pei_vasp_run_stretch.py``
     - 单轴 / 双轴应变扫描，拟合平衡晶格常数 ``a0``
   * - NEB
     - ``pei_vasp_run_neb``
     - 原子迁移能垒、缺陷扩散路径
   * - convergence
     - ``pei_vasp_run_convergence``
     - ENCUT / KPOINTS 收敛测试
   * - Cij
     - ``pei_vasp_run_cij_energy``
     - 二阶弹性常数（能量-应变法，5 个 Voigt 方向）
   * - cohesive
     - ``pei_vasp_run_cohesive``
     - cohesive energy、参考原子能量平台
   * - HOEC
     - ``pei_vasp_run_hoec_energy``
     - 高阶弹性常数（2/3/4 阶，Wang-Li 法）
   * - DOS / band
     - ``pei_vasp_run_dos_band``
     - 电子态密度与能带结构
   * - surface energy
     - ``pei_vasp_run_surface_energy``
     - 表面能（加真空层，bulk vs slab 对比）
   * - KPAR/NCORE
     - ``pei_vasp_run_kpar_ncore``
     - 并行效率基准（固定 workload，变 KPAR × NCORE）

输入文件结构
------------

所有 bulk workflow 统一从 ``y_full_relax`` 出发：它是一个已收敛的完整弛豫
结果，提供 ``CONTCAR`` （结构源）、``INCAR`` 、``KPOINTS`` 、``POTCAR`` 和
``sub*`` Slurm 脚本。每个 ``run_*`` 脚本在 ``y_full_relax`` 同级运行，把
生成的批量目录写到 ``y_<workflow>/y_dir/<case>`` 。

典型目录树（以 Au HCP 计算为例）：

.. code-block:: text

   project/
   |-- y_full_relax/                  # 统一入口：已收敛弛豫结果
   |   |-- CONTCAR                    # 结构源（POSCAR 由它复制）
   |   |-- INCAR                      # 弛豫参数（NSW=150, ISIF=3, IBRION=2）
   |   |-- KPOINTS                    # k 点网格（A 100 = Γ 中心 100×100×100）
   |   |-- POTCAR                     # 赝势（用户自备）
   |   `-- sub.544.sh                 # Slurm 提交脚本
   |-- y_stretch/                      # stretch workflow 输出
   |   |-- y_dir/
   |   |   |-- 0.99600000/            # 每个 case 一个目录
   |   |   |   |-- POSCAR INCAR KPOINTS POTCAR
   |   |   |   `-- sub.slurm
   |   |   |-- 0.99650000/
   |   |   `-- ... 17 个应变量
   |   |-- p_post_stretch.txt         # 拟合结果（a0, poly coeffs）
   |   `-- p_post_stretch.pdf         # 应变-能量曲线
   |-- y_cij_energy/
   |   |-- y_cij_energy_c11/y_dir/    # 5 个 Voigt 方向各一棵子树
   |   |-- y_cij_energy_c12/y_dir/
   |   |-- y_cij_energy_c13/y_dir/
   |   |-- y_cij_energy_c33/y_dir/
   |   |-- y_cij_energy_c44/y_dir/
   |   `-- y_post_cij_energy.txt      # 拟合 Cij + 派生量（E, ν, μ）
   `-- y_hoec_energy/
       |-- y_hoec_energy_M01/y_dir/   # hex 20 modes / cubic 11 modes
       |-- ...
       `-- y_post_hoec.json           # 2/3/4 阶常数（机器可读）

关键 INCAR 参数（来自 Au HCP ``y_full_relax`` ，仅展示参数部分）：

.. code-block:: bash

   ENCUT=   550        # 截断能 eV
   EDIFF=   1E-8       # 电子步收敛标准
   NSW=150             # 最大离子步
   ISIF= 3             # 3=全弛豫（cell+ions）; 2=只弛豫 ions
   IBRION= 2           # CG 法离子弛豫
   EDIFFG= -0.0005     # 力收敛标准 eV/Å
   ISMEAR= 1           # Methfessel-Paxton
   SIGMA=  0.2         # 展宽 eV
   NCORE= 4  NPAR=8    # 并行
   KPAR= 4             # k 点并行

提交前路径
----------

1. 用 ``mymetal.build``/ASE 生成结构，并检查元素、原子数、cell 和 PBC。
2. 在 ``y_src`` 准备 ``INCAR`` 、``KPOINTS`` 、``POTCAR`` 和结构源文件。
3. 运行 workflow 生成 ``y_dir/<case>`` ，逐项检查输入 diff。
4. 先让 Slurm 入口只生成脚本，运行 ``bash -n`` 并检查资源、module 和命令。
5. 明确添加 ``-if_sbatch`` 后才提交。
6. 分别核对 Slurm 状态、VASP 退出、电子/离子收敛，再做后处理。

安全检查
--------

.. code-block:: bash

   python slurm_utils/slurm_universal/pei_slurm_univ_submit.py -h
   bash -n vasp_utils/vasp_universal/pei_vasp_univ_sbatch
   bash -n vasp_utils/vasp_workflow_bulk/pei_vasp_run_neb

清理、续算和 resubmit 脚本可能删除输出或改写 POSCAR/INCAR。帮助文本不能替代
源码审查；先复制一个小型非生产目录验证。POTCAR 的来源、授权和元素顺序由用户
负责，本仓库文档不分发 POTCAR。

执行流程
--------

VASP workflow 的典型执行流程如下：

.. code-block:: text

   y_full_relax (已收敛)
        │
        ▼
   pei_vasp_run_<workflow>           ← 仅准备 y_<wf>/y_dir/<case>/
        │                               (cohesive 例外：自动 sbatch)
        ▼
   pei_slurm_univ_submit.py           ← 生成 sub_slurm_univ.sh
   -mode each_subdir -chunks K        (默认不提交，加 -if_sbatch 才提交)
        │
        ▼
   pei_vasp_univ_sbatch <job_dir>    ← 单目录 runner：srun vasp_std
   (exit 0=完成 / 1=失败 / 10=跳过)
        │
        ▼
   pei_vasp_univ_post                ← 批量后处理：刮 OUTCAR
        │                               → y_post_*.txt
        ▼
   pei_vasp_plot_all -<option>       ← 绘图：find y_<wf> → plot 脚本
                                        → p_post_*.pdf

Advanced workflow 概览
----------------------

* EOS：``vasp_utils/vasp_workflow_bulk/pei_vasp_run_eos.py`` 。
* NEB：``vasp_utils/vasp_workflow_bulk/pei_vasp_run_neb`` 与
  ``vasp_utils/neb_utils/`` 。
* GSFE：``vasp_utils/vasp_workflow_planar_defects/pei_vasp_run_gsfe`` 。
* decohesion：
  ``vasp_utils/vasp_workflow_planar_defects/pei_vasp_run_decohesion`` 。

这些入口真实存在，但第一轮没有经过许可和物理审核的最小 fixture，因此这里只
提供定位，不把它们包装成可复现 Tutorial。

实际结果展示（真实数据）
------------------------

以下图表均基于集群 (zcm6) 上 Au HCP 的真实 VASP 计算结果生成，数据来自
``construct_dataset/calculate/`` 下的 ``A11-1`` 、``A11-2`` 、``A11-3`` 和
``decohesion`` 子目录。

Stretch workflow 结果
~~~~~~~~~~~~~~~~~~~~~

面内双轴 (xy) 拉伸扫描，17 个应变点 (factor 0.996..1.004)，每点静态单点计算。
二次拟合得到平衡 factor 和 poly coefficients。

.. figure:: /_static/images/generated/vasp_stretch_real.png
   :alt: HCP Au xy 拉伸扫描真实结果
   :width: 90%

   Au HCP 面内 (xy) 拉伸扫描，17 个 VASP 静态单点（真实数据）。二次拟合
   平衡 factor = 0.99984，对应 a₀ = 2.85925 Å，c₀ = 4.80006 Å。

Cij workflow 结果
~~~~~~~~~~~~~~~~~

二阶弹性常数（能量-应变法），5 个 Voigt 方向 (C11, C12, C13, C33, C44)，
每方向 ~117 个应变点。拟合后得到 5 个独立常数及派生量 (E, ν, μ)。

.. figure:: /_static/images/generated/vasp_cij_real.png
   :alt: HCP Au Cij 真实结果
   :width: 95%

   Au HCP 二阶弹性常数（真实数据）。(a) C11 能量-应变拟合示意；
   (b) 5 个 Cij 的柱状图。C44 = 20.5 GPa 偏小，反映 HCP Au 的低剪切模量。

拟合结果汇总：

.. list-table::
   :header-rows: 1
   :widths: 25 25 50

   * - 常数
     - 值 (GPa)
     - 说明
   * - C11
     - 259.12
     - 面内压缩模量
   * - C12
     - 181.82
     - 面内耦合
   * - C13
     - 142.92
     - 面内-面外耦合
   * - C33
     - 242.75
     - c 方向压缩模量
   * - C44
     - 20.50
     - 基面剪切模量（偏小）
   * - E_x
     - 120.45
     - x 方向 Young's modulus
   * - E_z
     - 150.10
     - z 方向 Young's modulus
   * - ν_xy
     - 0.558
     - 面内 Poisson 比
   * - ν_xz
     - 0.260
     - 面内-面外 Poisson 比

HOEC workflow 结果
~~~~~~~~~~~~~~~~~~

高阶弹性常数（Wang-Li 能量-应变法），hex 对称性 20 个 mode，2/3/4 阶共 34
个独立常数。四阶对 ENCUT/KPOINTS/EDIFF 极其敏感，``y_full_relax`` 必须充分
收敛。

.. figure:: /_static/images/generated/vasp_hoec_real.png
   :alt: HCP Au 高阶弹性常数真实结果
   :width: 95%

   Au HCP 高阶弹性常数（真实数据）。二阶 5 个、三阶 10 个、四阶 19 个独立
   常数。四阶常数可达 ±77,000 GPa 量级，需要高精度参考计算。

Cohesive energy 结果
~~~~~~~~~~~~~~~~~~~~

cohesive energy 通过 scale 扫描 (0.60..4.00，42 个静态单点) 获得。当
scale=4.00 时原子几乎自由，该平台即参考原子能量 ``Eatom`` 。

.. figure:: /_static/images/generated/vasp_cohesive_real.png
   :alt: HCP Au cohesive energy 真实结果
   :width: 95%

   Au HCP cohesive energy（真实数据）。E₀ = -3.9169 eV/atom，
   E_atom = -0.0498 eV，E_coh = -3.8671 eV/atom，体弹模量 B = 188.17 GPa。

Decohesion 结果
~~~~~~~~~~~~~~~

decohesion（分离功）通过沿基面法向逐层分离扫描 (d=0..8.0 Å，39 个点) 获得。
γ(d) 曲线在 d ≈ 4.5 Å 处达到平台，平台值即理想分离功 γ∞。

.. figure:: /_static/images/generated/vasp_decohesion_real.png
   :alt: HCP Au decohesion 真实结果
   :width: 80%

   Au HCP decohesion 曲线（真实数据）。平台 γ∞ ≈ 1969 mJ/m²，对应理想
   解理能。

Stretch 收敛性
~~~~~~~~~~~~~~

拉伸扫描的平衡晶格常数随参与拟合的应变点数收敛：

.. figure:: /_static/images/generated/vasp_convergence_real.png
   :alt: 拉伸扫描收敛性
   :width: 75%

   拉伸扫描收敛性：随应变点数从 3 增到 17，拟合 a₀ 收敛到 server 报告值
   2.85925 Å。

参数说明
--------

Stretch workflow 参数（``pei_vasp_run_stretch.py`` ）：

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

HOEC workflow 参数（``pei_vasp_run_hoec_energy`` ）：

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - 参数
     - 说明
   * - ``-symmetry``
     - ``cubic`` 或 ``hex`` ，默认自动检测
   * - ``-emax``
     - 单轴参考模式（cubic A、hex M01）的最大振幅
   * - ``-de``
     - 单轴参考模式的步长
   * - ``-static``
     - ``NSW=0, IBRION=-1`` ，clamped-ion 常数
   * - ``-small_shear``
     - hex only：原位缩放剪切分量，保持亚稳 HCP branch
   * - ``-no_scale_window``
     - 所有 mode 共用一个窗口（旧行为，会过应变混合 mode）

与 mymetal 库函数的关联
-----------------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - 模块
     - 功能
   * - ``mymetal.io.vasp``
     - 读写 POSCAR/CONTCAR（workflow 脚本读取结构源）
   * - ``mymetal.post``
     - 解析 ``OUTCAR`` ，提取能量、力、应力、收敛状态
   * - ``mymetal.post.stretch``
     - stretch 后处理：拟合 poly coeffs，输出 a0
   * - ``mymetal.post.Cij_energy``
     - Cij 后处理：5 方向应变-能量拟合
   * - ``mymetal.post.gsfe``
     - GSFE 后处理：tilted-cell 剪切偏移扫描
   * - ``mymetal.calculate.calqm.kpoints``
     - k 点距离 → k 点网格（``get_size_by_distance`` ）

常见问题与排查
--------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - 问题
     - 排查
   * - OUTCAR 未收敛
     - 检查 ``OSZICAR`` 的 electronic/ionic 步数；用
       ``pei_vasp_univ_post`` 区分完成、未收敛、缺失与运行中
   * - stretch a0 偏移
     - 确认 ``-init/-final`` 对称覆盖平衡点；负应变列表用 ``-strains=...``
   * - C44 异常
     - C44 对应变窗口和拟合阶数敏感；检查 ``-e`` 参数和 poly 拟合 rms
   * - HOEC 四阶常数发散
     - ``y_full_relax`` 必须密 k 点 + 高 ENCUT + 紧 EDIFF (~1e-8)；
       检查 ``y_post_hoec_conv.pdf`` 的 plateau
   * - cohesive 无平台
     - 确认 scale 上限 ≥ 4.00；原子在 k=4 时应接近自由
   * - Slurm 作业超时
     - 用 ``-chunks K`` 拆分并发流；cohesive 的 42 点跑在一个 sequential
       job 里，wall time 要留足

后处理
------

少量目录可直接使用 :doc:`../tutorials/outcar_batch`。规模更大时，先用
``pei_vasp_univ_post`` 区分完成、未收敛、缺失与运行中，再将确认完成的
``OUTCAR`` 交给 :mod:`mymetal.post`。surface energy 的单位和等价表面数见
:doc:`../tutorials/surface_energy`。

.. seealso::

   - :doc:`vasp_universal` — 单目录 runner/post/清理工具
   - :doc:`vasp_workflow_bulk` — bulk 目录型 workflow（EOS、stretch、NEB、
     convergence、Cij、cohesive、DOS/band、surface energy、HOEC、KPAR/NCORE）
   - :doc:`neb_utils` — NEB 后处理工具集（``nebbarrier`` 、``nebef`` 、
     ``nebmovie`` 、``neb_plot`` ）
