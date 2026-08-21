文档维护
========

内容分工
--------

* ``getting_started/``：一条无需外部程序、可复制运行的新手路径。
* ``tutorials/``：围绕真实任务组织目标、输入、代码、输出、解释和检查。
* ``manual/``：解释多组件 workflow、状态与恢复边界。
* ``api/``：精确记录 Python 签名、单位、返回值和异常。
* ``reference/``：脚本、依赖与维护约定。

叙事说明使用中文；API 名、命令、文件名和原始输出保留英文。不要把未经执行的
片段写成已验证教程。

本地构建
--------

.. code-block:: bash

   python -m pip install -r docs/requirements.txt
   python -m pip install --no-deps "git+https://github.com/ShengLin1001/myalloy_package.git@master"
   python docs/scripts/generate_structure_images.py
   python -m sphinx -b html -W --keep-going docs/source docs/_build/html

完整最小验收
------------

.. code-block:: bash

   python docs/examples/getting_started_au111.py --output docs/_build/example-au111
   python docs/examples/surface_energy.py
   python docs/examples/outcar_summary.py mymetal/example/post-outcar-parse/y_dir \
       -cases 0.997 1.000 -output docs/_build/outcar-summary.csv
   python -m compileall mymetal docs/examples docs/scripts
   python -m sphinx -b html -W --keep-going docs/source docs/_build/html
   python -m sphinx -b linkcheck -W --keep-going docs/source docs/_build/linkcheck

Fixture 与图片
--------------

优先复用仓库内已有、小而稳定的 tracked fixture，并在教程中记录相对路径和
SHA-256。不要复制来源或许可不明确的数据。生成图必须由
``docs/scripts/`` 中的确定性脚本产生，并保留文字替代说明。

构建产物
--------

新构建统一写入被忽略的 ``docs/_build/``。历史 ``docs/build/`` 当前仍是 tracked
产物；文档升级不修改、覆盖或批量清理该目录。
