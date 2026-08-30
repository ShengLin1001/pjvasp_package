.. _tutorial-miller-index-and-density:

HCP Miller 指数转换与密度计算
==============================

:Audience: 需要转换 HCP 方向指数或计算结构密度的用户
:Time: 5 分钟
:Requires: Python、ASE、numpy、pjvasp_package
:Runs VASP: No

目标
----

演示两个 universal 工具：

1. :func:`mymetal.universal.atom.miller.three_index_to_four_index` —
   HCP 方向在三指数 ``[U,V,W]`` 和四指数 ``[u,v,t,w]`` 之间转换；
2. :func:`mymetal.universal.atom.density.cal_density` — 计算原子数密度
   ``natoms / volume``。

背景：HCP Miller 指数
----------------------

HCP 晶体方向常用四指数 ``[u,v,t,w]`` 表示，其中 ``t = -(u+v)``，以满足
三轴等价性。三指数 ``[U,V,W]`` 与四指数的转换关系：

.. math::

   u = \frac{2U - V}{3}, \quad v = \frac{2V - U}{3}, \quad t = -(u+v), \quad w = W

运行完整脚本
------------

.. code-block:: console

   $ python docs/examples/miller_index_and_density.py --output docs/_build/example-miller-density

.. literalinclude:: ../../examples/miller_index_and_density.py
   :language: python
   :linenos:

预期输出
--------

.. code-block:: text

   structure    formula  natoms  volume(A^3)   n_dens(at/A^3)  mass_dens(g/cm^3)
   ----------------------------------------------------------------------------
   FCC Cu       Cu4           4      47.0459         0.085023             8.9717
   BCC Fe       Fe2           2      23.6399         0.084603             7.8454
   HCP Mg       Mg2           2      46.4920         0.043018             1.7362
   Diamond Si   Si8           8     160.1030         0.049968             2.3303

结果图
------

.. figure:: /_static/images/generated/miller_index_and_density.png
   :alt: HCP Miller index conversion table and density bar chart for FCC/BCC/HCP/diamond

   左图：HCP 方向 3↔4 指数转换示例表。右图：FCC Cu、BCC Fe、HCP Mg、
   Diamond Si 的原子数密度和质量密度对比。

结果含义
--------

* **Miller 指数**： ``[1,0,0]`` → 四指数 ``[2/3, -1/3, -1/3, 0]`` ；
  ``[2,-1,0,0]`` → 三指数 ``[3,0,0]`` （与 ``[1,0,0]`` 同方向，
  round-trip 在整数缩放下一致）。
* **密度**：FCC Cu ``8.97 g/cm³`` 、BCC Fe ``7.85 g/cm³`` 与教科书值一致。
  HCP Mg 因 c/a≈1.62 且原子间距大，密度最低 ``1.74 g/cm³``。

验证方法
--------

* 3→4→3 round-trip 一致（整数缩放内）；
* 所有密度 > 0；
* 图片非空白。

相关 API
--------

* :func:`mymetal.universal.atom.miller.three_index_to_four_index`
* :func:`mymetal.universal.atom.density.cal_density`
* :doc:`neighbor_distances` （近邻距离分析）
