LAMMPS Workflows
================

本页是 LAMMPS workflow 的图文总览：从应用场景、输入文件结构、执行流程到基于
集群真实计算结果的展示。模板脚本的逐文件说明见 :doc:`lmp_utils`。

.. contents:: 本页内容
   :local:
   :depth: 2

应用场景
--------

``lmp_utils`` 是 LAMMPS 分子动力学计算的模板和 workflow 集合，用于材料力学
性质计算。当前仓库用它们作为以下 workflow 的起点：

.. list-table::
   :header-rows: 1
   :widths: 22 28 50

   * - workflow
     - 输出目录
     - 应用场景
   * - stretch
     - ``y_stretch/<lat>/``
     - 拉伸扫描，拟合平衡晶格常数 ``a0`` （FCC/BCC/HCP 三种晶格）
   * - Cij energy
     - ``y_Cij_energy/<lat>/``
     - 二阶弹性常数（能量-应变法，5 个 Voigt 方向）
   * - GSFE
     - ``y_gsfe/<lat>/<type>/``
     - 广义层错能（FCC: 111/100 面；HCP: 基面/柱面/锥面）

模板文件使用占位符（如 ``aa_template`` 、``pair_style_template`` ），由 runner
脚本在运行时通过 ``sed`` 替换为实际参数。模板与 ``mymetal`` Python package
刻意分离：生产运行前需根据目标势函数、unit style 和集群调度器进行配置。

.. note::

   这些模板是计算起点，不是即用型脚本。使用前需根据目标势函数、unit style
   和集群调度器进行配置。

输入文件结构
------------

LAMMPS workflow 的目录结构如下（以 Au UNNEP 势为例）：

.. code-block:: text

   project/
   |-- template/                     # LAMMPS 模板文件（含占位符）
   |   |-- stretch_template.in       # 拉伸计算入口
   |   |-- Cij_energy_template.in    # 弹性常数入口
   |   |-- gsfe_template.in          # GSFE 入口
   |   |-- general_init.mod          # units metal, boundary
   |   |-- general_potential.mod     # pair_style/pair_coeff（占位符）
   |   |-- general_mass.mod          # 原子质量（运行时写入）
   |   |-- stretch.mod              # 拉伸循环逻辑
   |   |-- stretch_full_relax.mod   # 完全弛豫
   |   |-- stretch_constrained_relax.mod  # 约束弛豫
   |   |-- Cij_energy.mod           # Cij 变形循环
   |   |-- gsfe.mod                 # GSFE 循环
   |   `-- gsfe_model.py            # GSFE 超胞构建（mymetal.build）
   |-- potential/                    # 势函数文件
   |   `-- UNEP-v1-main.txt          # UNNEP 势（用户自备）
   |-- y_stretch/                    # stretch 输出
   |   |-- fcc/
   |   |   |-- y_dir/                # 101 个应变量目录
   |   |   |-- dump/                 # 每步 dump 文件
   |   |   |-- p_post_stretch.txt    # 拟合结果（a0, poly coeffs）
   |   |   `-- p_post_stretch.pdf    # 应变-能量曲线
   |   |-- bcc/
   |   `-- hcp/
   |-- y_Cij_energy/
   |   |-- fcc/
   |   |   |-- y_dir/
   |   |   |-- dump/
   |   |   |-- y_post_cij_energy.txt  # 拟合 Cij + 派生量
   |   |   `-- y_post_cij_energy.pdf
   |   |-- bcc/
   |   `-- hcp/
   |-- y_gsfe/
   |   |-- fcc/
   |   |   |-- FCC_111/
   |   |   |   |-- y_dir/
   |   |   |   |-- dump/
   |   |   |   |-- y_post_gsfe.txt
   |   |   |   `-- y_post_gsfe.pdf
   |   |   `-- FCC_100/
   |   `-- hcp/
   |       |-- HCP_basal/
   |       |-- HCP_prism1w/
   |       |-- HCP_pyr1w/
   |       `-- HCP_pyr2/
   |-- post/                         # 后处理脚本
   |   |-- stretch.py
   |   |-- Cij_energy.py
   |   |-- gsfe.py
   |   `-- post.ipynb
   `-- sub_slurm_univ.sh             # Slurm 提交脚本

关键模板内容
~~~~~~~~~~~~

``stretch_template.in`` （拉伸入口）：

.. code-block:: bash

   include general_init.mod          # units metal, boundary
   include stretch_model.mod          # 晶格参数 (aa, lat)
   include general_potential.mod      # pair_style/pair_coeff
   include stretch_full_relax.mod     # 完全弛豫
   include stretch.mod                # 拉伸循环

``stretch.mod`` 中的拉伸循环核心（101 步，±0.3%）：

.. code-block:: bash

   variable step equal 101
   variable large_strain equal 0.003
   variable i loop ${step}
   variable id equal ${i}-(${step}-1)/2-1    # -50 .. +50
   variable factor equal ${id}/((${step}-1)/2)  # -1 .. +1
   variable strain equal ${large_strain}*${factor}
   variable stretch equal ${strain}+1
   change_box all x delta 0 ${deltaxx} y delta 0 ${deltayy} z delta 0 ${deltazz}

执行流程
--------

LAMMPS workflow 的典型执行流程：

.. code-block:: text

   template/ (含占位符) + potential/ (势函数)
        │
        ▼
   pei_lmp_run_properties [pair_style] [pair_coeff] [python_path] [mass] [template]
        │  (sed 替换占位符 → 生成实际 in.* 文件)
        │
        ├─ stretch (fcc/bcc/hcp) → y_stretch/<lat>/ → post/stretch.py
        │   ↓ 读取平衡晶格常数 a0
        ├─ Cij_energy (fcc/bcc/hcp) → y_Cij_energy/<lat>/ → post/Cij_energy.py
        │   ↓ 同上
        ├─ gsfe (fcc: 111/100, hcp: basal/prism1w/pyr1w/pyr2)
        │   → y_gsfe/<lat>/<type>/ → post/gsfe.py
        │
        └─ Summary (汇总三个 post 的退出码)

runner 带有 ``srun`` step 创建失败的自动重试机制（``LMP_MAX_TRY`` 次，默认 99
次，每次间隔 ``LMP_RETRY_SLEEP`` 秒，默认 5s），只重试瞬时调度错误。

初始晶格常数由 runner 内置：

.. list-table::
   :header-rows: 1
   :widths: 20 30 50

   * - 晶格
     - latnum
     - 初始 a (Å)
   * - hcp
     - 1
     - 2.88
   * - fcc
     - 2
     - 4.08
   * - bcc
     - 3
     - 3.20

实际结果展示（真实数据）
------------------------

以下图表均基于集群 (zcm6) 上 Au UNNEP 势的真实 LAMMPS 计算结果，数据来自
``pj-test-properties-gold/`` 目录。

Stretch workflow 结果
~~~~~~~~~~~~~~~~~~~~~

FCC Au 各向同性 (xyz) 拉伸，101 个应变点 (factor 0.997..1.003)，二次拟合
得到平衡晶格常数。

.. figure:: /_static/images/generated/lammps_stretch_real.png
   :alt: FCC Au 拉伸扫描真实结果
   :width: 90%

   FCC Au 各向同性 (xyz) 拉伸扫描（真实数据）。二次拟合平衡 factor =
   1.00001，对应 a₀ = 4.15852 Å，E = -3.22695 eV/atom。

Cij workflow 结果
~~~~~~~~~~~~~~~~~

FCC Au 二阶弹性常数（能量-应变法），5 个 Voigt 方向。由于 FCC 立方对称性，
C13 = C12，C33 = C11。

.. figure:: /_static/images/generated/lammps_cij_real.png
   :alt: FCC Au Cij 真实结果
   :width: 95%

   FCC Au 二阶弹性常数（真实数据）。(a) C11 能量-应变拟合示意；
   (b) 5 个 Cij 柱状图。C44 = 32.7 GPa。

拟合结果汇总：

.. list-table::
   :header-rows: 1
   :widths: 25 25 50

   * - 常数
     - 值 (GPa)
     - 说明
   * - C11
     - 154.72
     - 立方对称：= C33
   * - C12
     - 118.32
     - 立方对称：= C13
   * - C13
     - 118.32
     - 立方对称：= C12
   * - C33
     - 154.72
     - 立方对称：= C11
   * - C44
     - 32.69
     - 剪切模量
   * - E_x
     - 52.17
     - Young's modulus
   * - ν_xy
     - 0.433
     - Poisson 比
   * - μ_xz
     - 32.69
     - 剪切模量

GSFE workflow 结果
~~~~~~~~~~~~~~~~~~

FCC Au (111) 面广义层错能，21 个剪切偏移点，用 tilted-cell method 施加
剪切偏移。局部极大值对应不稳定层错能 (USFE)。

.. figure:: /_static/images/generated/lammps_gsfe_real.png
   :alt: FCC Au (111) GSFE 真实结果
   :width: 85%

   FCC Au (111) 广义层错能曲线（真实数据）。局部极大 γ ≈ 106.5 mJ/m²
   （不稳定层错能），Asf = 7.488 Å²。

GSFE 滑移系类型与剪切偏移参数：

.. list-table::
   :header-rows: 1
   :widths: 20 15 15 50

   * - gsfe_type
     - bp1
     - bp2
     - 说明
   * - ``FCC_111``
     - -0.5
     - 1/3
     - FCC (111) 面
   * - ``FCC_100``
     - 0.5
     - 0.5
     - FCC (100) 面
   * - ``HCP_basal``
     - -0.5
     - 1/3
     - HCP 基面
   * - ``HCP_prism1w``
     - 0.5
     - 0
     - HCP 柱面 I
   * - ``HCP_pyr1w``
     - 0
     - 0.5
     - HCP 锥面 I
   * - ``HCP_pyr2``
     - 0
     - 0.5
     - HCP 锥面 II

参数说明
--------

``pei_lmp_run_properties`` 参数（按位置）：

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   * - 参数
     - 默认值
     - 说明
   * - ``pair_style``
     - ``eam``
     - LAMMPS pair_style
   * - ``pair_coeff``
     - ``./potential/Au_u3.eam``
     - 势函数文件路径
   * - ``python_path``
     - ``python``
     - Python 可执行文件路径
   * - ``mass_content``
     - ``mass 1 196.97``
     - 原子质量（写入 ``general_mass.mod`` ）
   * - ``lmp_template_path``
     - ``./template``
     - 模板文件目录

sed 模板替换占位符：

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - 占位符
     - 替换为
     - 出现在
   * - ``aa_template``
     - 初始晶格常数 (fcc:4.08, bcc:3.20, hcp:2.88)
     - stretch_model.mod, gsfe_model.py
   * - ``cc_template``
     - HCP 的 c 参数
     - gsfe_model.py
   * - ``lat_template``
     - latnum (1=hcp, 2=fcc, 3=bcc)
     - stretch_model.mod, gsfe_model.py
   * - ``pair_style_template``
     - pair_style (如 eam)
     - general_potential.mod
   * - ``pair_coeff_template``
     - pair_coeff 文件路径
     - general_potential.mod
   * - ``gsfe_type_template``
     - 滑移系类型
     - gsfe_model.py, gsfe_template.in
   * - ``gsfe_bp1_template``
     - bp1 剪切偏移
     - gsfe.mod
   * - ``gsfe_bp2_template``
     - bp2 剪切偏移
     - gsfe.mod
   * - ``python_path_template``
     - Python 路径
     - gsfe_template.in

与 mymetal 库函数的关联
-----------------------

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - 模块
     - 功能
   * - ``mymetal.build.bulk.gsfe``
     - 构建 GSFE 超胞（``create_gsfe_model`` 写出 ``data.ini`` ）
   * - ``mymetal.post.stretch``
     - stretch 后处理（``post_lammps_stretch`` 拟合 a0）
   * - ``mymetal.post.Cij_energy``
     - Cij 后处理（``post_lammps_Cij_energy`` 拟合弹性常数）
   * - ``mymetal.post.gsfe``
     - GSFE 后处理（``post_gsfe`` 计算 γ(s)）

后处理脚本
----------

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - 脚本
     - 功能
   * - ``post/stretch.py``
     - 调用 ``mymetal.post.stretch.post_lammps_stretch`` ，生成
       ``p_post_stretch.txt`` 和 ``p_post_stretch.pdf``
   * - ``post/Cij_energy.py``
     - 调用 ``mymetal.post.Cij_energy.post_lammps_Cij_energy`` ，生成
       ``y_post_cij_energy.txt`` 和 ``y_post_cij_energy.pdf``
   * - ``post/gsfe.py``
     - 调用 ``mymetal.post.gsfe.post_gsfe`` ，生成 ``y_post_gsfe.txt`` 、
       ``y_post_gsfe.pdf`` 和 ``y_post_gsfe.u3.pdf``

常见问题与排查
--------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - 问题
     - 排查
   * - LAMMPS 启动失败
     - runner 自动重试 ``LMP_MAX_TRY`` 次（默认 99）；检查 module 和
       ``pair_coeff`` 路径
   * - stretch a0 偏移
     - 确认 ``stretch_full_relax.mod`` 完全弛豫后再扫描；检查
       ``large_strain`` 和 ``step``
   * - Cij 非对称
     - FCC 应满足 C13=C12、C33=C11；偏差大说明应变窗口或拟合有问题
   * - GSFE 曲线不光滑
     - 检查 ``gsfe_model.py`` 超胞原子数；tilted-cell method 对 cell
       精度敏感
   * - 势函数文件找不到
     - 确认 ``pair_coeff_template`` 占位符被正确替换；路径用绝对路径

.. seealso::

   - :doc:`lmp_utils` — LAMMPS 模板脚本逐文件说明
   - :doc:`vasp_workflow_bulk` — VASP 对应的 stretch/Cij/GSFE workflow
   - :doc:`slurm_utils` — Slurm 提交入口
   - :doc:`workflows` — Workflow 总览
   - :doc:`../tutorials/gsfe_models` — GSFE 模型构建教程
   - :doc:`../tutorials/cij_energy_fitting` — 弹性常数拟合教程
   - :doc:`../tutorials/biaxial_stretch` — 双轴应变教程
