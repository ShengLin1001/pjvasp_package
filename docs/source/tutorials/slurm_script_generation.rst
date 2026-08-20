.. _tutorial-slurm-script-generation:

Slurm 作业脚本生成（dry-run）
==============================

:Audience: 需要生成 Slurm 脚本但不立即提交的用户
:Time: 5 分钟
:Requires: Python、pjvasp_package
:Runs VASP: No
:Submits SLURM job: No

目标
----

演示 :mod:`mymetal.slurm.submit` 的脚本生成能力（dry-run 模式），
不调用 ``sbatch``。展示从 ``generate_script_header`` 到完整脚本组装的
流程。

.. warning::

   本教程只生成脚本文本，不提交任何 SLURM 作业。``if_sbatch=False``
   是默认行为。

运行完整脚本
------------

.. code-block:: console

   $ python docs/examples/slurm_script_generation.py --output docs/_build/example-slurm

.. literalinclude:: ../../examples/slurm_script_generation.py
   :language: python
   :linenos:

预期输出
--------

.. code-block:: text

   ===== generate_script_header =====
   #SBATCH --partition=debug
   #SBATCH --nodes=1
   #SBATCH --ntasks=4
   #SBATCH --time=01:00:00
   ...

   ===== generate_slurm_script_base =====
   #!/bin/bash
   #SBATCH --partition=debug
   ...
   pei_slurm_univ_launch_retry srun -n $SLURM_NTASKS  vasp_std

   ===== check_wall_time =====
   '01:00:00' -> valid
   '2-00:00:00' -> valid
   'invalid' -> caught SystemExit

   wrote: .../slurm_script_generation.png
   OK: all assertions passed.

结果图
------

.. figure:: /_static/images/generated/slurm_script_generation.png
   :alt: Slurm script generation pipeline flowchart and rendered script text

   左图：脚本生成流程图（header → launcher → assemble → write）。
   右图：生成的完整 Slurm 脚本文本渲染。

验证方法
--------

* 生成的脚本包含 ``#SBATCH`` 行；
* 包含 partition 名称；
* ``check_wall_time`` 对合法输入返回字符串，对非法输入抛 ``SystemExit``；
* 图片非空白。

相关 API
--------

* :func:`mymetal.slurm.submit.pei_slurm_univ_submit`
* :func:`mymetal.slurm.submit.check_wall_time`
* :doc:`../manual/slurm` （SLURM 工作流完整指南）
* :doc:`../reference/scripts` （pei_* 脚本参考）
