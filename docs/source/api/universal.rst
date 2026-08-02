Universal Utilities
-------------------

.. automodule:: mymetal.universal
   :members:

Atom operations
---------------

.. automodule:: mymetal.universal.atom.fixatom
   :members:
   :show-inheritance:

.. automodule:: mymetal.universal.atom.moveatom
   :members:
   :show-inheritance:

.. automodule:: mymetal.universal.atom.neighbor
   :members:
   :show-inheritance:

``three_index_to_four_index``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.atom.miller.three_index_to_four_index

用途
   HCP 晶体方向在三指数 ``[U,V,W]`` 和四指数 ``[u,v,t,w]``（``t=-(u+v)``）
   之间转换。``reverse=False`` 时 3→4，``reverse=True`` 时 4→3。

Related tutorials
   :doc:`../tutorials/miller_index_and_density`

``cal_density``
~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.atom.density.cal_density

用途
   返回原子数密度 ``natoms / volume``，单位 atoms/Å³。

Related tutorials
   :doc:`../tutorials/miller_index_and_density`

.. automodule:: mymetal.universal.atom.rotate
   :members:
   :show-inheritance:

.. automodule:: mymetal.universal.atom.delatom
   :members:
   :show-inheritance:

Checks and transformations
--------------------------

.. automodule:: mymetal.universal.check.checkinput
   :members:
   :show-inheritance:

.. automodule:: mymetal.universal.check.convergence
   :members:
   :show-inheritance:

.. automodule:: mymetal.universal.check.type
   :members:
   :show-inheritance:

``get_cna_count`` / ``check_phase_transition``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.check.atoms.get_cna_count

.. autofunction:: mymetal.universal.check.atoms.get_cna_count_dict

.. autofunction:: mymetal.universal.check.atoms.check_phase_transition

用途
   基于 OVITO Common Neighbor Analysis 统计 FCC/HCP/BCC/ICO/OTHER 相分布。
   ``check_phase_transition`` 比较初始和最终结构的 CNA 分布，判断是否发生
   相变。

.. note::

   这些函数依赖可选的 ``ovito`` 包。

.. automodule:: mymetal.universal.data.dataadjust
   :members:
   :show-inheritance:

.. automodule:: mymetal.universal.data.datachange
   :members:
   :show-inheritance:

.. automodule:: mymetal.universal.data.patterntrans
   :members:
   :show-inheritance:

.. automodule:: mymetal.universal.math.operations
   :members:
   :show-inheritance:

.. automodule:: mymetal.universal.matrix.adjust
   :members:
   :show-inheritance:

Plotting
--------

``periodic_table_heatmap`` / ``van_arkel_triangle``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.plotting.periodic_table_heatmap
   :no-index:

.. autofunction:: mymetal.universal.plot.plotting.van_arkel_triangle
   :no-index:

用途
   ``periodic_table_heatmap`` 在周期表上按元素值着色（如体弹模量、结合能）。
   ``van_arkel_triangle`` 把一组材料按电负性绘制在 van Arkel 三角图上，
   区分离子/共价/金属键特征。

Related tutorials
   :doc:`../tutorials/periodic_table_and_arkel`

.. note::

   ``mymetal.universal.plot.plotting`` 还包含 ``pretty_plot``、
   ``pretty_polyfit_plot``、``format_formula`` 等通用绘图工具。完整 API
   见源码；此处仅展示两个特色函数。

.. automodule:: mymetal.universal.plot.workflow
   :members:
   :show-inheritance:
   :exclude-members: my_plot_neb_xy, my_plot_neb, my_plot_neb_full
