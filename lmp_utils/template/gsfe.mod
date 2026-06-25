############################ gsfe loop ############################


################# Check odd step #################
# Custom set parameters, for cubic plot
# v_step must be odd number

variable step equal 21
# 每个滑移系的剪切方向/大小都不同（见 vasp_utils/.../pei_vasp_run_gsfe 的 bp1/bp2）。
# DFT 把 c 矢量的 x、y 分量分别平移 bp1*a1、bp2*b2；在 LAMMPS 中 c 矢量的 x、y 分量
# 就是 xz、yz tilt，因此 Δxz = bp1*lx、Δyz = bp2*ly。runner 按 gsfe_type 注入下面两个值。
variable bp1 equal gsfe_bp1_template
variable bp2 equal gsfe_bp2_template

variable check_even equal ${step}%2

if "${check_even} < 1e-10" then &
    "variable temp equal ${step}+1" &
else &
    "variable temp equal ${step}" &

variable step equal ${temp}

################## Check end ##################


################### loop ###################
label loop
clear
read_restart ../dump/restart.equil
include general_potential.mod

# Get stretch factor v_stretch
variable i  loop ${step}                  # 1, 2, ..., step
variable id equal ${i}                   # 1, 2, ..., step
variable factor  equal (${id}-1)/(${step}-1)  # 0, 1/10, 2/10, ..., 1.

# lammps cell
# xx 0  0
# xy yy 0
# xz yz zz
# 按 gsfe_type 的 bp1/bp2 同时剪 xz、yz，复刻 DFT 的偏移约定（与 pei_vasp_run_gsfe 一致）：
#   Δxz = bp1 * lx * factor,  Δyz = bp2 * ly * factor,  factor: 0 -> 1
# $(lx)/$(ly) 立即取平衡盒子（未倾斜）的边长，避免被本步 tilt 影响。
# 先把盒子升格为三斜：FCC_100 等正交超胞 ASE 写出无 tilt，LAMMPS 会拒绝直接加 tilt
# (ERROR: Cannot change box tilt factors for orthogonal box)；change_box triclinic 仅把
# domain->triclinic 置 1，对已是三斜的盒子(HCP_basal/FCC_111)无副作用。
# lx = xhi - xlo = xx, ly = yhi - ylo = yy
variable deltaxz equal ${bp1}*$(lx)*${factor}
variable deltayz equal ${bp2}*$(ly)*${factor}
change_box all triclinic
change_box all xz delta ${deltaxz} yz delta ${deltayz}  units box

# gsfe: \sigma33 = 0, full relax in z position, and fixed in x,y position
# red point: \sigma3j = 0, full relax in x, y, z position
# Here is gsfe setting
fix restrict  all setforce 0 0 NULL
fix br all box/relax z  0   fixedpoint 0 0 0 scaleyz no scalexy no scalexz no
min_style   cg
minimize ${etol} ${ftol} ${maxiter} ${maxeval}
unfix       br
unfix restrict


variable jobn  string "$(v_id-1:%16.12e)"
include general_output.mod

next i
jump SELF loop
