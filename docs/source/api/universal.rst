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

.. automodule:: mymetal.universal.search.find
   :members:
   :show-inheritance:

.. automodule:: mymetal.universal.print.print
   :members:
   :show-inheritance:

.. automodule:: mymetal.universal.print.printafter
   :members:
   :show-inheritance:

用途
   ``pr`` / ``er`` / ``warn`` / ``fail`` 是带 emoji 前缀的输出助手；
   ``fail`` 会 raise ``SystemExit(1)``。``find_line_position`` /
   ``extract_line_at_position`` 按标记从文本文件提取行，供 OUTCAR 等
   parser 使用。

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
   :doc:`../tutorials/periodic_table_and_arkel`、
   :doc:`../tutorials/plot_gallery_plotting`

.. note::

   ``mymetal.universal.plot.plotting`` 还包含 ``pretty_plot``、
   ``pretty_polyfit_plot``、``format_formula`` 等通用绘图工具。完整 API
   见源码；此处仅展示两个特色函数。

``my_plot`` / ``my_plot_colorbar`` / ``my_plot_modify_ploted_figure``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.plot.my_plot

.. autofunction:: mymetal.universal.plot.plot.my_plot_colorbar

.. autofunction:: mymetal.universal.plot.plot.my_plot_modify_ploted_figure

用途
   画布创建器：``my_plot`` 是标准入口，返回配好样式的 ``(fig, axes)``；
   ``my_plot_colorbar`` 创建带右侧 colorbar 坐标区的画布，返回
   ``(fig, (ax_main, ax_cbar))``；``my_plot_modify_ploted_figure`` 改造
   已画好的图（常用于 DOS）。

Related tutorials
   :doc:`../tutorials/plot_gallery_plot`、
   :doc:`../tutorials/plot_gallery`

``general_set_all_rcParams`` / ``general_modify_legend`` / 标注助手
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.general.general_set_all_rcParams

.. autofunction:: mymetal.universal.plot.general.general_modify_legend

.. autofunction:: mymetal.universal.plot.general.add_color_band

.. autofunction:: mymetal.universal.plot.general.add_circle_number

.. autofunction:: mymetal.universal.plot.general.add_arrow

.. autofunction:: mymetal.universal.plot.general.generate_gradient_colors

.. autofunction:: mymetal.universal.plot.general.general_adjust_text

用途
   通用样式与标注助手：``general_set_all_rcParams`` 一次性配置全部
   rcParams 并返回图例修饰闭包；``general_modify_legend`` 统一图例边框；
   ``add_color_band`` / ``add_circle_number`` / ``add_arrow`` 在图上叠加
   色带、编号圆圈、箭头；``generate_gradient_colors`` 生成渐变色；
   ``general_adjust_text`` 自动避让文字重叠。

Related tutorials
   :doc:`../tutorials/plot_gallery_general`、
   :doc:`../tutorials/plot_gallery`

``my_plot_convergence`` / ``my_plot_relax_convergence`` / ``my_plot_stretch``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.workflow.my_plot_convergence
   :no-index:

.. autofunction:: mymetal.universal.plot.workflow.my_plot_relax_convergence
   :no-index:

.. autofunction:: mymetal.universal.plot.workflow.my_plot_stretch
   :no-index:

.. autofunction:: mymetal.universal.plot.workflow.my_plot_kpar_ncore
   :no-index:

用途
   工作流后处理出图：收敛测试、离子弛豫、拉伸拟合、KPAR/NCORE 基准。

Related tutorials
   :doc:`../tutorials/plot_gallery_workflow`、
   :doc:`../tutorials/plot_gallery`

``my_plot_energy_components``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.energy.my_plot_energy_components

用途
   绘制 OUTCAR 能量分量随自变量的变化并做多项式拟合。

Related tutorials
   :doc:`../tutorials/plot_gallery_energy`、
   :doc:`../tutorials/plot_gallery`

``my_plot_interlayer_distance`` / ``my_plot_zpositions`` / ``my_plot_rdf``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.atominfo.my_plot_interlayer_distance

.. autofunction:: mymetal.universal.plot.atominfo.my_plot_zpositions

.. autofunction:: mymetal.universal.plot.atominfo.my_plot_rdf

用途
   原子结构信息绘图：层间距、z 位置、径向分布函数。

Related tutorials
   :doc:`../tutorials/plot_gallery_atominfo`、
   :doc:`../tutorials/plot_gallery`

``my_plot_learning_curve`` / ``my_plot_compare`` / ``my_plot_rmse_by_tag``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.n2p2.my_plot_learning_curve

.. autofunction:: mymetal.universal.plot.n2p2.my_plot_compare

.. autofunction:: mymetal.universal.plot.n2p2.my_plot_rmse_by_tag

.. autofunction:: mymetal.universal.plot.n2p2.my_plot_epoch_rmse

用途
   n2p2 训练诊断：学习曲线、DFT-vs-NNP 散点、按 tag 的 RMSE 柱状图、
   epoch 监控。

Related tutorials
   :doc:`../tutorials/plot_gallery_n2p2`、
   :doc:`../tutorials/plot_gallery`

``my_render``
~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.render.my_render

用途
   用 OVITO 管线渲染结构图片。

Related tutorials
   :doc:`../tutorials/plot_gallery_render`、
   :doc:`../tutorials/plot_gallery`

``ppt2picture``
~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.ppt.ppt2picture

用途
   将 PowerPoint 文件的全部幻灯片导出为裁剪后的图片。

Related tutorials
   :doc:`../tutorials/plot_gallery_ppt`、
   :doc:`../tutorials/plot_gallery`

DOS and band plotting
~~~~~~~~~~~~~~~~~~~~~

.. automodule:: mymetal.universal.plot.oldplotdos
   :members:
   :show-inheritance:

用途
   旧版 DOS/能带绘图工具：``my_plot_complete_dos`` 绘制总态密度，
   ``my_plot_idos`` 绘制分波态密度，``my_plot_element_spd_dos`` 绘制
   元素 s/p/d 分轨道 DOS。依赖可选的 ``pymatgen`` DOS 数据结构。

Related tutorials
   :doc:`../tutorials/plot_gallery_oldplotdos`、
   :doc:`../tutorials/plot_gallery`

.. automodule:: mymetal.universal.plot.workflow
   :members:
   :show-inheritance:
   :exclude-members: my_plot_neb_xy, my_plot_neb, my_plot_neb_full
