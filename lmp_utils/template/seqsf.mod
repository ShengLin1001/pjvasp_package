############################ seqsf: faulted stacking sequence ############################
# 换结构只能靠 clear 重来。clear 会把 units/atom_style/boundary 打回默认，
# 所以 read_data 之前必须重新 include general_init.mod；变量本身不受 clear 影响，
# 但重新 include 会把 etol/ftol 等 equal-style 变量按同值覆盖，无副作用。
clear
include general_init.mod
read_data ./data.fault
include general_potential.mod

# 理想层错构型（弛豫前），供 post 算 Yin 式 (1) 的原子位移 u = X - x
write_dump all custom ../dump/ideal.lammpstrj id type x y z modify sort id format line "%2d %2d %16.12e %16.12e %16.12e"

# 与参考胞同一套设置：原子全放开，σ3j = 0。对称性会让 xz/yz 停在 0，
# 放开它们只是为了和 ssf.mod 保持同一口径，不改变结果。
fix br all box/relax z  0  xz  0  yz  0  fixedpoint 0 0 0 scaleyz no scalexy no scalexz no
min_style   cg
minimize ${etol} ${ftol} ${maxiter} ${maxeval}
unfix       br

variable i    equal 2
variable jobn string "1.000000000000e+00"
include general_output.mod
