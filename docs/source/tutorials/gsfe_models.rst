.. _tutorial-gsfe-models:

构建 FCC/HCP 广义层错能（GSFE）模型
=====================================

:Audience: 需要为 GSFE 计算准备超胞结构的用户
:Time: 5 分钟
:Requires: Python、ASE、pjvasp_package
:Runs VASP: No

目标与边界
----------

演示 :func:`mymetal.build.bulk.gsfe.create_gsfe_model` 如何根据
``gsfe_type`` 一次性生成面向 GSFE 计算的超胞。本教程构建三种代表性取向：
FCC(111)、HCP(0001) 基面和 HCP(10-10) prism I（wide），打印结构摘要并渲染
侧视图对比图。不运行 VASP，不调用 ``sbatch``。

.. note::

   ``create_gsfe_model`` 支持六种 ``gsfe_type``。其中 ``"FCC_111"``、
   ``"HCP_basal"``、``"HCP_prism1w"`` 仅依赖 ASE；``"HCP_pyr1w"``、
   ``"HCP_pyr2"``、``"FCC_100"`` 路径调用 ``vasp_create_*`` 系列函数，
   依赖可选的 ``myvasp`` 包。本教程只展示纯 ASE 路径。

背景
----

广义层错能（GSFE）描述晶体沿特定滑移面发生刚性剪切位移时的能量变化，
是理解位错核心结构和滑移系活性的关键量。不同晶体结构的典型滑移面：

* FCC：(111) 面，``a/2<110>`` 全位错，``a/6<112>`` Shockley 不全位错；
* HCP：(0001) 基面、``{10-10}`` prism I、``{10-11}`` pyramid I/II。

``create_gsfe_model`` 为每种取向生成对应超胞，后续通过在面内施加位移并
做 VASP 静态计算即可扫描 GSFE 曲线。

运行完整脚本
------------

.. code-block:: console

   $ python docs/examples/gsfe_models.py --output docs/_build/example-gsfe

.. literalinclude:: ../../examples/gsfe_models.py
   :language: python
   :linenos:

预期输出
--------

.. code-block:: text

   label gsfe_type formula atoms a_A b_A c_A alpha beta gamma volume_A3 pbc
   Au FCC(111) FCC_111 Au18 18 2.8850 2.8850 42.4006 90.00 90.00 120.00 305.6279 [True, True, True]
   Mg HCP(0001) basal HCP_basal Mg12 12 3.2100 3.2100 31.2600 90.00 90.00 120.00 278.9521 [True, True, True]
   Mg HCP(10-10) prism I (wide) HCP_prism1w Mg24 24 3.2100 5.2100 33.3593 90.00 90.00 90.00 557.9042 [True, True, True]
   wrote: .../docs/_build/example-gsfe/gsfe_models.png

结果图
------

.. figure:: /_static/images/generated/gsfe_models.png
   :alt: Side views of FCC(111), HCP(0001) basal and HCP(10-10) prism I GSFE supercells

   三个 GSFE 超胞的侧视图。左：FCC Au(111)，7 层，18 原子；
   中：HCP Mg(0001) 基面，10 层，12 原子；右：HCP Mg(10-10) prism I
   （wide），5 层，24 原子。三方向 PBC 均为 True。

结果含义
--------

* FCC(111) 超胞 ``c ≈ 42.4 Å``，包含 7 层 Au，每层 3 个原子（primitive
  面内单元），层间距 ``a/sqrt(3) ≈ 2.36 Å``。
* HCP(0001) 基面超胞 ``c ≈ 31.3 Å``，包含 10 层 Mg，层间距 ``c/2 ≈ 2.6 Å``。
* HCP(10-10) prism I（wide）超胞 ``b ≈ 5.21 Å``、``c ≈ 33.4 Å``，
  wide 层间距 ``a/sqrt(3) ≈ 1.85 Å``，共 24 个原子。

验证方法
--------

脚本在写出结果前执行以下断言：

* 三种模型的原子数均大于 0；
* 化学式分别为 ``Au18``、``Mg12``、``Mg24``；
* PBC 三方向均为 ``True``；
* 图片文件存在且非近似纯白。

常见错误
--------

* ``ImportError: myvasp``：使用了 ``"HCP_pyr1w"`` / ``"HCP_pyr2"`` /
  ``"FCC_100"`` 类型但未安装 ``myvasp`` 包。改用纯 ASE 类型，或安装
  companion package（见 :doc:`../user_guide/install`）。
* ``a`` 或 ``c`` 为 ``None``：HCP 类型必须同时提供 ``a`` 和 ``c``。

下一步
------

* 用 :func:`mymetal.build.film.stretch.stretch_along_direction_to_cell`
  对 GSFE 超胞施加面内位移扫描，见 :doc:`cubic_cell_and_stretch`。
* GSFE 后处理（拟合曲线、提取 USF/USFE）见 :doc:`../manual/workflows`
  中的 Advanced Workflow 概览。

相关 API
--------

* :func:`mymetal.build.bulk.gsfe.create_gsfe_model`
* :func:`mymetal.build.bulk.create.create_fcc_111`
* :func:`mymetal.build.bulk.create.create_hcp_basal`
* :func:`mymetal.build.bulk.create.create_hcp_prism1`
