lmp_utils LAMMPS 工作流
========================

``lmp_utils`` 是 LAMMPS 分子动力学计算的模板和 workflow 集合，用于材料力学
性质计算：拉伸 (stretch)、弹性常数 (Cij energy) 和广义层错能 (GSFE)。
模板文件使用占位符，由 runner 脚本在运行时通过 ``sed`` 替换为实际参数。

.. note::

   这些模板是计算起点，不是即用型脚本。使用前需根据目标势函数、unit style
   和集群调度器进行配置。

.. contents:: 本页内容
   :local:
   :depth: 2

整体架构
--------

``pei_lmp_run_properties`` 是主 runner，依次执行三种计算：

.. code-block:: text

   pei_lmp_run_properties [pair_style] [pair_coeff] [python_path] [mass_content] [lmp_template_path]
       │
       ├─ stretch (fcc/bcc/hcp)  → y_stretch/<lat>/  → post/stretch.py
       │   ↓ 读取平衡晶格常数 a0, c0
       ├─ Cij_energy (fcc/bcc/hcp) → y_Cij_energy/<lat>/ → post/Cij_energy.py
       │   ↓ 同上
       ├─ gsfe (fcc: FCC_111/FCC_100, hcp: HCP_basal/prism1w/pyr1w/pyr2)
       │   → y_gsfe/<lat>/<type>/ → post/gsfe.py
       │
       └─ Summary (汇总三个 post 的退出码)

runner 带有 ``srun`` step 创建失败的自动重试机制（``LMP_MAX_TRY`` 次，
默认 99 次，每次间隔 ``LMP_RETRY_SLEEP`` 秒，默认 5s），只重试瞬时调度错误。

lmp_universal
-------------

pei_lmp_run_properties
~~~~~~~~~~~~~~~~~~~~~~

主 runner 脚本，在一个 Slurm 分配内执行完整的 stretch → Cij → gsfe 流程。

.. code-block:: bash

   # 默认：Au EAM 势
   pei_lmp_run_properties

   # 自定义势函数和质量
   pei_lmp_run_properties eam /path/to/Au_u3.eam python "mass 1 196.97" /path/to/template

参数（按位置）：

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
     - 原子质量 (写入 ``general_mass.mod``)
   * - ``lmp_template_path``
     - ``./template``
     - 模板文件目录

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

template 模板文件
-----------------

三种计算各有独立的 ``.in`` 入口模板，通过 ``include`` 引入共用的 ``.mod`` 模块。

stretch_template.in
~~~~~~~~~~~~~~~~~~~

拉伸计算模板。完全弛豫后施加面内/面外拉伸扫描。

.. code-block:: bash

   include general_init.mod          # units metal, boundary
   include stretch_model.mod          # 晶格参数 (aa, lat)
   include general_potential.mod      # pair_style/pair_coeff
   include stretch_full_relax.mod     # 完全弛豫
   include stretch.mod                # 拉伸循环

Cij_energy_template.in
~~~~~~~~~~~~~~~~~~~~~~~

弹性常数 energy-strain 法模板。约束弛豫后对 5 个 Voigt 方向施加变形。

.. code-block:: bash

   include general_init.mod
   include stretch_model.mod
   include general_potential.mod
   include stretch_constrained_relax.mod   # 约束弛豫
   # 5 个变形方向: C11, C12, C13, C33, C44
   variable dirn equal 1
   include Cij_energy.mod
   # ... (dirn = 2..5)

gsfe_template.in
~~~~~~~~~~~~~~~~~

广义层错能模板。完全弛豫后用 tilted-cell method 施加剪切偏移。

.. code-block:: bash

   include general_init.mod
   shell python gsfe_model.py        # 生成 data.ini
   read_data ./data.ini
   include general_potential.mod
   minimize ...                      # 弛豫
   include gsfe.mod                  # GSFE 循环 (change_box)

.mod 模块文件
~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - 文件
     - 功能
   * - ``general_init.mod``
     - ``units metal``, ``boundary p p p``, ``atom_style atomic``
   * - ``general_mass.mod``
     - 原子质量 （运行时由 runner 写入）
   * - ``general_potential.mod``
     - ``pair_style`` / ``pair_coeff`` （运行时 ``sed`` 替换占位符）
   * - ``general_output.mod``
     - 输出格式定义
   * - ``general_structural_info.mod``
     - cell 信息提取
   * - ``stretch_model.mod``
     - 晶格参数变量 (``aa_template``, ``lat_template``)
   * - ``stretch.mod``
     - 拉伸循环逻辑
   * - ``stretch_full_relax.mod``
     - 完全弛豫 (fcc/bcc: jump bcc_fcc, hcp: jump hcp)
   * - ``stretch_constrained_relax.mod``
     - 约束弛豫 (固定 cell 形状)
   * - ``gsfe.mod``
     - GSFE 循环 (``gsfe_bp1_template``, ``gsfe_bp2_template``)
   * - ``Cij_energy.mod``
     - Cij 变形循环
   * - ``Cij_energy_output.mod``
     - Cij 输出格式

gsfe_model.py
~~~~~~~~~~~~~

用 ``mymetal.build.bulk.gsfe.create_gsfe_model`` 构建 GSFE 超胞，写出
``data.ini`` (LAMMPS data 格式)。模板中的占位符由 runner ``sed`` 替换：

- ``aa_template`` → 晶格常数 a
- ``cc_template`` → 晶格常数 c
- ``lat_template`` → latnum (1=hcp, 2=fcc, 3=bcc)
- ``gsfe_type_template`` → 滑移系类型

GSFE 滑移系类型与剪切偏移：

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

post 后处理脚本
---------------

post/stretch.py
~~~~~~~~~~~~~~~

调用 ``mymetal.post.stretch.post_lammps_stretch``，处理
``y_stretch/<lat>/dump/`` 下的数据，生成 ``p_post_stretch.txt`` 和
``p_post_stretch.pdf``。

对 fcc/bcc/hcp 三种晶格循环处理。

post/Cij_energy.py
~~~~~~~~~~~~~~~~~~~

调用 ``mymetal.post.Cij_energy.post_lammps_Cij_energy``，处理
``y_Cij_energy/<lat>/dump/`` 下的数据，拟合弹性常数，生成
``y_post_cij_energy.txt`` 和 ``y_post_cij_energy.pdf``。

post/gsfe.py
~~~~~~~~~~~~

调用 ``mymetal.post.gsfe.post_gsfe``，处理
``y_gsfe/<lat>/<type>/dump/`` 下的数据，生成 ``y_post_gsfe.txt``、
``y_post_gsfe.pdf`` 和 ``y_post_gsfe.u3.pdf``。

支持的滑移系：

- FCC: ``FCC_111``, ``FCC_100``
- HCP: ``HCP_basal``, ``HCP_prism1w``, ``HCP_pyr1w``, ``HCP_pyr2``

sed 模板替换
------------

``pei_lmp_run_properties`` 在运行时用 ``sed`` 替换模板中的占位符：

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

示例脚本
--------

.. literalinclude:: ../../examples/lmp_utils_overview.py
   :language: python
   :caption: lmp_utils 子包功能总览示例

.. seealso::

   - :doc:`lammps` — LAMMPS 简版说明
   - :doc:`vasp_workflow_bulk` — VASP 对应的 stretch/Cij/GSFE workflow
   - :doc:`slurm_utils` — Slurm 提交入口
   - :doc:`workflows` — Workflow 总览
   - :doc:`../tutorials/gsfe_models` — GSFE 模型构建教程
   - :doc:`../tutorials/cij_energy_fitting` — 弹性常数拟合教程
   - :doc:`../tutorials/biaxial_stretch` — 双轴应变教程
