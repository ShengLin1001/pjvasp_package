.. _tutorial-plot-gallery-energy:

energy 模块 —— 能量分量分解
============================

:模块: :mod:`mymetal.universal.plot.energy`
:函数数: 1

模块概述
--------

``energy`` 模块用于绘制 VASP OUTCAR 中的能量分量（TOTEN、PSCENC、TEWEN、
DENC、EXHF、XCENC 等）随某自变量（如应变、目录索引）的变化，并对每个分量
做多项式拟合。

示例图
------

.. figure:: /_static/images/generated/plot_gallery_energy.png
   :width: 960px
   :alt: my_plot_energy_components demonstration

   用合成 OUTCAR 风格字典绘制的能量分量面板，每个子图对应一个能量 tag，
   数据点 + 多项式拟合曲线。

公开函数
--------

my_plot_energy_components
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.energy.my_plot_energy_components
   :no-index:

用途
   生成能量分量面板图。输入 ``dic_list`` 是字典列表，每个字典必须包含
   ``directory`` 键和全部 OUTCAR 能量 tag 键（小写）。

.. note::

   函数内部遍历固定 tag 列表（来自 ``tagsref``），包括 ``exhf``、
   ``xcenc``、``double1``、``double2``、``double``、``ediel_sol`` 等。如果
   字典缺少任何一个 tag，会抛出 ``KeyError``。

最小示例
   .. code-block:: python

      from mymetal.universal.plot.energy import my_plot_energy_components

      dirs = [0.0, 0.1, 0.2, 0.3]
      tags = ["toten", "pscenc", "tewen", "denc", "exhf", "xcenc",
              "double1", "double2", "double", "eentro", "ebands",
              "eatom", "ediel_sol"]
      dic = {"directory": dirs}
      for t in tags:
          dic[t] = [base[t] * (1 + 0.05 * x) for x in dirs]
      my_plot_energy_components(
          dic_list=[dic], label_list=["set A"], atoms_number=4,
          save_path="energy.png", poly_fit_e=2,
      )

相关 API
--------

* :mod:`mymetal.universal.plot.energy`
* :doc:`plot_gallery_workflow` （工作流出图）
* :doc:`plot_gallery` （总览）
