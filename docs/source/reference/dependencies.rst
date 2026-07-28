依赖边界
========

Python runtime
--------------

常用路径依赖 ``numpy``、``scipy``、``ase``、``spglib``、``pandas`` 和
``matplotlib``。完整版本声明以仓库根目录 ``requirements.txt`` 和
``setup.py`` 为准；文档构建使用较小的 ``docs/requirements.txt``。

Companion package
-----------------

``myalloy`` 和可导入的 ``myvasp`` 来自
``https://github.com/ShengLin1001/myalloy_package``：

.. code-block:: bash

   python -m pip install "git+https://github.com/ShengLin1001/myalloy_package.git@master"

本仓库的 ``myvasp/`` 目录主要保存 shell/VTST 脚本，不是上述 Python
package 的替代品。

可选 Python 依赖
----------------

``hetbuilder``、``ovito``、``pymatgen``、``torch`` 和 ``torchvision``
只服务于部分模块。仅在调用对应功能时安装；不要为 Au(111) Getting Started
安装整套可选依赖。

外部程序与 HPC 服务
-------------------

VASP、LAMMPS、n2p2、Slurm、MPI launcher、环境 module 和 POTCAR 均不是
Python 依赖，也不会由 ``pip`` 安装。目标集群上先检查：

.. code-block:: bash

   module avail
   command -v sbatch
   command -v srun

仓库脚本
--------

``vasp_utils/``、``slurm_utils/``、``lmp_utils/``、``n2p2_utils/`` 和
``myvasp/`` 是源码树的一部分。项目尚未把它们注册为 ``console_scripts``；
请使用明确路径或配置 ``PATH``，详见 :doc:`scripts`。
