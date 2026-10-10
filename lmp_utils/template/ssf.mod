############################ ssf: single tilted point ############################
# 层错矢量 (t1, t2)，单位为面内基矢 a1/a2 的分数，runner 按 <gsfe_type>_<fault> 注入。
# 与 gsfe.mod 完全一致的映射：DFT 把 c 矢量平移 (t1*a1, t2*b2)，LAMMPS 里就是
#   Δxz = t1 * lx,  Δyz = t2 * ly
# 初猜只需对称性给出的近似位置：xz/yz 是自由度，弛豫会把层错矢量推到真实的 t + rn。
variable t1 equal ssf_t1_template
variable t2 equal ssf_t2_template

# $(lx)/$(ly) 立即求值，取的是参考胞（未倾斜）的边长，避免被本步 tilt 影响。
variable deltaxz equal ${t1}*$(lx)
variable deltayz equal ${t2}*$(ly)

# 正交盒子不能直接加 tilt（ERROR: Cannot change box tilt factors for orthogonal box），
# 先升格三斜；change_box triclinic 对已是三斜的盒子（HCP_basal / pyr1）无副作用。
change_box all triclinic
change_box all xz delta ${deltaxz} yz delta ${deltayz}  units box

# 理想层错构型（弛豫前）。post 侧按 Yin 式 (1) 用它算 u = X - x，
# 并据此判断这个层错到底有没有面内弛豫（有 => 标准法会算错）。
write_dump all custom ../dump/ideal.lammpstrj id type x y z modify sort id format line "%2d %2d %16.12e %16.12e %16.12e"

# SSF 的定义性设置，与 gsfe.mod 的两处差别：
#   - 没有 fix setforce：原子 x, y, z 全放开；
#   - box/relax 放开 z + xz + yz：σ33 = σ13 = σ23 = 0。
# 面内 (N1a1, N2a2) 仍然锁死在体相值 —— Yin §3 强调，否则层错的本征应力被胞内其余
# 部分的反向应力抵消，弹性能随 N3 漂移，算出来的层错能没有意义。
fix br all box/relax z  0  xz  0  yz  0  fixedpoint 0 0 0 scaleyz no scalexy no scalexz no
min_style   cg
minimize ${etol} ${ftol} ${maxiter} ${maxeval}
unfix       br

variable i    equal 2
variable jobn string "1.000000000000e+00"
include general_output.mod
