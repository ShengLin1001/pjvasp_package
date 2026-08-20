.. _tutorial-plot-gallery-ppt:

ppt 模块 —— PowerPoint 幻灯片导出
=================================

:模块: :mod:`mymetal.universal.plot.ppt`
:函数数: 1

模块概述
--------

``ppt`` 模块通过 Windows PowerPoint COM 接口将 ``.ppt`` / ``.pptm`` /
``.pptx`` 文件的每一张幻灯片导出为裁剪后的图片文件。支持 PNG / JPG /
BMP / GIF / TIFF 等格式，自动裁剪幻灯片背景白边并保留指定内边距。

依赖：Microsoft PowerPoint + ``pywin32`` （仅 Windows）。Pillow 用于图片
裁剪和格式验证。

示例图
------

.. figure:: /_static/images/generated/plot_gallery_ppt.png
   :width: 720px
   :alt: ppt2picture call signature illustration

   ``ppt2picture`` 的调用签名和导出流程示意。

公开函数
--------

ppt2picture
~~~~~~~~~~~

.. autofunction:: mymetal.universal.plot.ppt.ppt2picture
   :no-index:

用途
   将一个 PowerPoint 文件的全部幻灯片导出为裁剪后的图片文件。

最小示例
   .. code-block:: python

      from mymetal.universal.plot.ppt import ppt2picture

      paths = ppt2picture(
          path_ppt="slides.pptx",
          path_output="./slides_images",
          style="png", dpi=600, padding_mm=1.0,
      )

相关 API
--------

* :mod:`mymetal.universal.plot.ppt`
* :doc:`plot_gallery` （总览）
