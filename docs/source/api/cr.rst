Contact-Resonance (CR) Analysis
================================

.. automodule:: mymetal.cr

Stiffness fitting (crsum3)
--------------------------

.. note::

   ``mymetal.cr.crsum3`` 从 CR（接触共振）测量数据拟合每条曲线的刚度。
   通过 Hertz 风格幂律曲线拟合 force-distance 数据，计算 CR 频率。
   该模块在文档构建时无法自动导入（依赖运行时环境）；
   函数签名和参数详见源码 ``mymetal/cr/crsum3.py``。

k-contact graph plotting (crplotkcontactgraph)
----------------------------------------------

.. note::

   ``mymetal.cr.crplotkcontactgraph`` 读取 CR 测量参数，绘制 k-contact
   graph（刚度接触图），用于可视化不同接触点的刚度分布。
   该模块在文档构建时无法自动导入；函数签名详见源码
   ``mymetal/cr/crplotkcontactgraph.py``。

.. note::

   CR 分析模块依赖实验数据文件，当前不提供无需数据的可运行示例。
   函数签名和参数详见源码或 :doc:`../manual/workflows` 中的 CR 概览。
