安装与环境
==========

建议使用 Python 3.10，并在独立虚拟环境中安装。以下命令只安装 Python
package；VASP、LAMMPS、n2p2、Slurm 和 POTCAR 均需由目标 HPC 环境另行提供。

安装 Python package
-------------------

.. code-block:: bash

   git clone https://github.com/ShengLin1001/pjvasp_package.git
   cd pjvasp_package
   python -m pip install -r requirements.txt
   python -m pip install -e .

验证导入和版本：

.. code-block:: bash

   python -c "import mymetal; print(mymetal.__file__)"
   python -c "import ase, matplotlib, numpy; print(ase.__version__)"

Companion package
-----------------

部分历史 VASP helper 会导入 ``myalloy`` 或 ``myvasp``。它们由 companion
repository 提供，而不是本仓库中的 ``myvasp/`` 脚本目录：

.. code-block:: bash

   python -m pip install "git+https://github.com/ShengLin1001/myalloy_package.git@master"

脚本入口与 PATH
---------------

当前 ``setup.py`` 没有注册 ``console_scripts``。仓库脚本应使用明确路径运行，
或把对应目录加入 ``PATH``：

.. code-block:: bash

   python slurm_utils/slurm_universal/pei_slurm_univ_submit.py -h
   bash vasp_utils/vasp_workflow_bulk/pei_vasp_run_neb -h

   export PATH="$PWD/slurm_utils/slurm_universal:$PWD/vasp_utils/vasp_universal:$PATH"
   command -v pei_slurm_univ_submit.py

CentOS HPC 上若命令缺失或版本不合适，先查看集群 module：

.. code-block:: bash

   module avail
   command -v python
   command -v sbatch

最小功能验收
------------

Getting Started 不依赖 VASP、POTCAR 或 Slurm：

.. code-block:: bash

   python docs/examples/getting_started_au111.py --output docs/_build/example-au111

成功时会写出 ``POSCAR`` 和 ``au111_slab.png``，并完成 ASE round-trip 检查。

构建文档
--------

.. code-block:: bash

   python -m pip install -r docs/requirements.txt
   python -m pip install --no-deps "git+https://github.com/ShengLin1001/myalloy_package.git@master"
   python -m sphinx -b html -W --keep-going docs/source docs/_build/html

生成文件统一放在 ``docs/_build/``；仓库中历史 ``docs/build/`` 不参与当前构建。
