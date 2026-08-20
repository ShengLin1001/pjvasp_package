.. _tutorial-plot-gallery-render:

render 模块 —— OVITO 管线渲染
==============================

:模块: :mod:`mymetal.universal.plot.render`
:函数数: 1

模块概述
--------

``render`` 模块封装 OVITO 的渲染管线，将 ``Pipeline`` 对象用指定相机方向
和渲染器输出为图片。依赖可选的 ``ovito`` 包。

.. note::

   OVITO 渲染需要真实的结构文件和 OVITO 安装，超出本 VASP-free 画廊的
   范围。下图是调用签名示意面板。

示例图
------

.. figure:: /_static/images/generated/plot_gallery_render.png
   :width: 720px
   :alt: my_render call signature illustration

   ``my_render`` 的调用签名和典型渲染流程示意。

公开函数
--------

my_render
~~~~~~~~~

.. autofunction:: mymetal.universal.plot.render.my_render
   :no-index:

用途
   用指定相机方向和渲染器渲染 OVITO 管线为图片。

最小示例
   .. code-block:: python

      from ovito.io import import_file
      from ovito.modifiers import CommonNeighborAnalysisModifier
      from ovito.vis import TachyonRenderer
      from mymetal.universal.plot.render import my_render

      pipeline = import_file('POSCAR')
      pipeline.modifiers.append(CommonNeighborAnalysisModifier())
      my_render(
          pipeline=pipeline, imagefile='render.png',
          size=(3200, 2400),
          renderer=TachyonRenderer(shadows=True, direct_light_intensity=1.1),
          camera_dir=(1, 1, 1),
      )

相关 API
--------

* :mod:`mymetal.universal.plot.render`
* :doc:`plot_gallery` （总览）
