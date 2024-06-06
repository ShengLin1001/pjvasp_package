http://bbs.keinsci.com/thread-19985-1-1.html

extforce.sh脚本的思路很简单：1）从OUTCAR中提取所有的TOTAL-FORCE部分信息（每个原子上的force vectors），
2）用for循环调用screenforce.py筛选出每一个离子步受力大于EDIFFG的原子，及对应的原子序号，存储在lt_ediffg_${i}.dat，
i为相应的离子步数。然后可以用p4vasp软件查看 具体是哪些原子没有达到力收敛标准。

使用的时候需要在extforce.sh里面设置一下自己体系的原子数（N_atoms），在screenforce.py里设置一下自己的EDIFFG（ediffg）。
需用python3，如果只有python2的话可以把screenforce.py的print那几行改成py2的格式即可。

使用时将extforce.sh, screenforce.py, OUTCAR放在同一目录下，执行./extforce.sh即可。