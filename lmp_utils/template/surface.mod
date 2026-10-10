############################ surface: cleave + relax ############################
# 真空厚度必须 > 2*rcut，否则两个表面隔着周期边界还能互相看见。
# 本项目 HDNNP 的截断半径是 6.01 Ang，取 20 Ang 真空有充分余量；runner 按需注入。
variable vacuum equal surface_vacuum_template

# 刚性劈开：只把 zhi 抬高，xz/yz/xy 保持不变（面内基矢锁死在体相值）。
# change_box z delta 的两个值是 (zlo 增量, zhi 增量)。
change_box all z delta 0 ${vacuum} units box

# 未弛豫表面：单点，不动原子
run 0
variable i    equal 2
variable jobn string "1.000000000000e+00"
include general_output.mod

# 弛豫表面：原子全放开，胞完全锁死。
# 这里不能加 box/relax：zz 已经含真空，放开它只会把真空压回去。
min_style   cg
minimize ${etol} ${ftol} ${maxiter} ${maxeval}
variable i    equal 3
variable jobn string "2.000000000000e+00"
include general_output.mod
