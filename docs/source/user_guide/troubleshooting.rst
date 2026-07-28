故障排查
========

``myvasp`` 或 ``myalloy`` 无法导入
----------------------------------

导入模块由 companion repository 提供：

.. code-block:: bash

   python -m pip install "git+https://github.com/ShengLin1001/myalloy_package.git@master"
   python -c "import myalloy, myvasp; print(myalloy.__file__)"

本仓库的 ``myvasp/`` 是 workflow 脚本目录，不等价于 importable package。

``pei_*: command not found``
----------------------------

本项目未注册 console script。先用仓库内路径确认脚本存在，再配置 ``PATH``：

.. code-block:: bash

   python slurm_utils/slurm_universal/pei_slurm_univ_submit.py -h
   export PATH="$PWD/slurm_utils/slurm_universal:$PATH"
   command -v pei_slurm_univ_submit.py

在 CentOS HPC 上，安装替代软件前先执行 ``module avail``。

POSCAR/CONTCAR 读回不一致
-------------------------

确认输入文件坐标模式、元素顺序和原子数；随后做最小 round-trip：

.. code-block:: python

   from mymetal.io.vasp import my_read_vasp, my_write_vasp

   atoms, scale = my_read_vasp("CONTCAR")
   my_write_vasp("POSCAR.check", atoms, lattice_scale_factor=scale)
   print(len(atoms), atoms.get_chemical_formula(), atoms.cell)

``OUTCAR`` 汇总为空
-------------------

``PostData`` 需要真实 ``OUTCAR``。先确认路径和关键行：

.. code-block:: bash

   test -s y_dir/1.000/OUTCAR
   grep -m 1 "Elapsed time" y_dir/1.000/OUTCAR

未收敛、文件截断和路径错误是不同状态，不应统一当作零值。

Slurm 生成、排队与计算失败
--------------------------

* 没有 ``-if_sbatch``：只生成脚本，这是预期 dry-run。
* ``PENDING``：检查 partition、account、资源和排队原因。
* ``RUNNING`` 后失败：检查作业 stdout、launcher、module 和程序退出码。
* 父作业失败：先确认子作业是否仍在队列，避免重复提交。

文档构建警告
------------

本手册将 warning 视为错误。先运行：

.. code-block:: bash

   python -m sphinx -b html -W --keep-going docs/source docs/_build/html

若 autodoc 只缺少不参与文档执行的可选依赖，可在 ``conf.py`` 的
``autodoc_mock_imports`` 中记录；真实教程依赖不能用 mock 掩盖。
