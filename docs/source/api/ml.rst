Machine Learning Utilities
--------------------------

.. automodule:: mymetal.ml
   :members:

n2p2 data and symmetry functions
--------------------------------

.. automodule:: mymetal.ml.n2p2.dataset
   :members:
   :show-inheritance:

.. automodule:: mymetal.ml.n2p2.calculate.sf

.. autofunction:: mymetal.ml.n2p2.calculate.sf.generate_radial_blocks

.. autofunction:: mymetal.ml.n2p2.calculate.sf.generate_angular_blocks

用途
   生成 n2p2 ACSF（atom-centered symmetry functions）的 radial 和 angular
   symmetry function 参数块。``generate_radial_blocks`` 生成径向 SF，
   ``generate_angular_blocks`` 生成角 SF。

.. automodule:: mymetal.ml.n2p2.sfparamgen.sfparamgen

CUR 特征选择
~~~~~~~~~~~~

.. autofunction:: mymetal.ml.n2p2.calculate.cur.cur_select

用途
   用 CUR 分解从 symmetry function 特征矩阵中选择最有代表性的子集，
   降低 n2p2 训练的特征冗余。

.. automodule:: mymetal.ml.n2p2.calculate.post
   :members:
   :show-inheritance:
   :exclude-members: build_rmse_by_tag_df

n2p2 training workflow
~~~~~~~~~~~~~~~~~~~~~~

.. note::

   ``mymetal.ml.n2p2.workflow`` 提供 ``PeiN2p2`` 类，封装 n2p2 的 dataset
   生成、symmetry function 配置、训练和后处理完整流程。它依赖 n2p2 可执行
   文件和 ``myvasp`` 包。该模块的 docstring 当前包含 RST 格式问题，暂不
   自动渲染完整 API；类和方法签名详见源码 ``mymetal/ml/n2p2/workflow.py``。

Optional image-model helpers
----------------------------

.. automodule:: mymetal.ml.dataset
   :members:
   :show-inheritance:

.. automodule:: mymetal.ml.model
   :members:
   :show-inheritance:

.. note::

   ``mymetal.ml.dataset`` 和 ``mymetal.ml.model`` 依赖 ``torch`` /
   ``torchvision`` （可选），提供图像分类模型的训练和可视化工具。

.. automodule:: mymetal.ml.confusionmatrix
   :members:
   :show-inheritance:

.. automodule:: mymetal.ml.plot
   :members:
   :show-inheritance:
