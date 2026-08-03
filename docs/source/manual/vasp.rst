VASP workflow
=============

功能分层
--------

``mymetal.io.vasp``
   Python 层读写 POSCAR/CONTCAR；不运行 VASP。

``mymetal.post``
   解析已存在的 ``OUTCAR`` 与批量目录；不判断调度器状态。

``vasp_utils/vasp_universal``
   单目录 runner、monitor、post、续算和清理脚本。

``vasp_utils/vasp_workflow_*``
   EOS、stretch、NEB、surface energy、GSFE 等目录型 workflow。

``myvasp``
   历史 shell 与 VTST-style 工具。可导入的 ``myvasp`` Python helper 来自
   ``myalloy_package``。

提交前路径
----------

1. 用 ``mymetal.build``/ASE 生成结构，并检查元素、原子数、cell 和 PBC。
2. 在 ``y_src`` 准备 ``INCAR``、``KPOINTS``、``POTCAR`` 和结构源文件。
3. 运行 workflow 生成 ``y_dir/<case>``，逐项检查输入 diff。
4. 先让 Slurm 入口只生成脚本，运行 ``bash -n`` 并检查资源、module 和命令。
5. 明确添加 ``-if_sbatch`` 后才提交。
6. 分别核对 Slurm 状态、VASP 退出、电子/离子收敛，再做后处理。

安全检查
--------

.. code-block:: bash

   python slurm_utils/slurm_universal/pei_slurm_univ_submit.py -h
   bash -n vasp_utils/vasp_universal/pei_vasp_univ_sbatch
   bash -n vasp_utils/vasp_workflow_bulk/pei_vasp_run_neb

清理、续算和 resubmit 脚本可能删除输出或改写 POSCAR/INCAR。帮助文本不能替代
源码审查；先复制一个小型非生产目录验证。POTCAR 的来源、授权和元素顺序由用户
负责，本仓库文档不分发 POTCAR。

Advanced workflow 概览
----------------------

* EOS：``vasp_utils/vasp_workflow_bulk/pei_vasp_run_eos.py``。
* NEB：``vasp_utils/vasp_workflow_bulk/pei_vasp_run_neb`` 与
  ``vasp_utils/neb_utils/``。
* GSFE：``vasp_utils/vasp_workflow_planar_defects/pei_vasp_run_gsfe``。
* decohesion：
  ``vasp_utils/vasp_workflow_planar_defects/pei_vasp_run_decohesion``。

这些入口真实存在，但第一轮没有经过许可和物理审核的最小 fixture，因此这里只
提供定位，不把它们包装成可复现 Tutorial。

后处理
------

少量目录可直接使用 :doc:`../tutorials/outcar_batch`。规模更大时，先用
``pei_vasp_univ_post`` 区分完成、未收敛、缺失与运行中，再将确认完成的
``OUTCAR`` 交给 :mod:`mymetal.post`。surface energy 的单位和等价表面数见
:doc:`../tutorials/surface_energy`。

``vasp_universal`` 的完整脚本说明见 :doc:`vasp_universal`。
