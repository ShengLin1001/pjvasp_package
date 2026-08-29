# Lattice parameters
# initial lattice constant
variable aa equal aa_template
variable lat equal lat_template   # 1 - "hcp", 2 - "fcc"(cubic), 3 - "bcc", 4 - "fcc_hex"
variable shift equal 1e-5

# ${lat_hex} = 1 表示六方框架（a2 与 a1 成 120°、盒子带 xy 倾斜、c 轴独立）。
# 下游模板（md_heating.mod 的 box/relax 与 npt）靠它决定用 "couple xy + z 独立"
# 还是 "iso"，而不是去硬编码某个 lat 编号 —— 否则每加一个六方相都要改一遍下游。
variable lat_hex equal 0

# ${aa} 的含义随 lat 变：立方相（2/3）是**常规立方胞的晶格常数**，
# 六方相（1/4）是**面内最近邻距离**。两者差 √2 倍，传参时别弄混。
print "Lattice type: ${lat}"

if "${lat} == 1" then &
    "jump SELF hcp"     &
elif "${lat} == 2"    &
    "jump SELF fcc"     &
elif "${lat} == 3"    &
    "jump SELF bcc"     &
elif "${lat} == 4"    &
    "jump SELF fcc_hex" &
else &
    "print 'Unknown lattice type: ${lat}'"

# BCC
label bcc
variable a11 equal  1.0
variable a22 equal  1.0
variable a33 equal  1.0

variable a21 equal  0.0
variable a31 equal  0.0
variable a32 equal  0.0

lattice custom ${aa}   &
                 a1    ${a11}              0.00000000000000       0.00000000000000 &
                 a2    0.00000000000000    ${a22}                 0.00000000000000 &
                 a3    0.00000000000000    0.00000000000000       ${a33}           &
                 basis 0.00000000000000    0.00000000000000       0.00000000000000 &
                 basis 0.50000000000000    0.50000000000000       0.50000000000000

jump SELF end_if

# FCC
label fcc
variable a11 equal  1.0
variable a22 equal  1.0
variable a33 equal  1.0

variable a21 equal  0.0
variable a31 equal  0.0
variable a32 equal  0.0

lattice custom ${aa}   &
                 a1    ${a11}              0.00000000000000       0.00000000000000 &
                 a2    0.00000000000000    ${a22}                 0.00000000000000 &
                 a3    0.00000000000000    0.00000000000000       ${a33}           &
                 basis 0.00000000000000    0.00000000000000       0.00000000000000 &
                 basis 0.50000000000000    0.50000000000000       0.00000000000000 &
                 basis 0.50000000000000    0.00000000000000       0.50000000000000 &
                 basis 0.00000000000000    0.50000000000000       0.50000000000000

jump SELF end_if

# HCP
label hcp
variable lat_hex delete
variable lat_hex equal 1
variable a11 equal  1.0
variable a22 equal 0.86602540378444
variable a33 equal 1.63299316185545

variable a21 equal -0.5
variable a31 equal 0.0
variable a32 equal 0.0

lattice custom ${aa}   &
                 a1    ${a11}              0.00000000000000       0.00000000000000 &
                 a2    ${a21}              ${a22}                 0.00000000000000 &
                 a3    0.00000000000000    0.00000000000000       ${a33}           &
                 basis 0.00000000000000    0.00000000000000       0.00000000000000 &
                 basis 0.33333333333333    0.66666666666666       0.50000000000000

jump SELF end_if

# FCC，但以密排面 (111) 为底的六方胞（3 原子 / 胞，ABC 堆垛）
#
# 为什么要这一支：常规立方 FCC 胞是 4 原子，最小超胞 2x2x2 = 32 原子；
# 换成 (111) 为底的 3 原子胞后 2x2x2 = 24 原子，与 HCP 的 2x2x3 = 24 完全一致，
# 也与初始训练集的原子数口径（[20,30] 逼近 24）对齐，DFT 代价直接降三成。
# 附带好处：FCC 与 HCP 从此共用同一个六方框架、同一个面内 ${aa}，
# 两相的盒子形状与 k 点网格都可比，相稳定性的对照更干净。
#
# 几何：${aa} = 面内最近邻距离（= a_cubic/sqrt(2)）。
#   a3/aa = sqrt(6) = 2.449...：3 个 (111) 层的总厚度 = 3 * a_cubic/sqrt(3) = sqrt(3)*a_cubic，
#                               除以 aa = a_cubic/sqrt(2) 即 sqrt(6)。
#   basis 的三个横向偏移 (0,0) -> (1/3,2/3) -> (2/3,1/3) 就是 A -> B -> C 三个不同的
#   三角子格位（与 HCP 的 A -> B -> A 只差第三层落在哪里）。
label fcc_hex
variable lat_hex delete
variable lat_hex equal 1
variable a11 equal  1.0
variable a22 equal 0.86602540378444
variable a33 equal 2.44948974968170

variable a21 equal -0.5
variable a31 equal 0.0
variable a32 equal 0.0

lattice custom ${aa}   &
                 a1    ${a11}              0.00000000000000       0.00000000000000 &
                 a2    ${a21}              ${a22}                 0.00000000000000 &
                 a3    0.00000000000000    0.00000000000000       ${a33}           &
                 basis 0.00000000000000    0.00000000000000       0.00000000000000 &
                 basis 0.33333333333333    0.66666666666666       0.33333333333333 &
                 basis 0.66666666666666    0.33333333333333       0.66666666666666

jump SELF end_if

label end_if
variable a11_box equal ${aa}*${a11}
variable a22_box equal ${aa}*${a22}
variable a33_box equal ${aa}*${a33}

variable a21_box equal ${aa}*${a21}  # xy
variable a31_box equal ${aa}*${a31}  # xz
variable a32_box equal ${aa}*${a32}  # yz

# Origin: (xlo,ylo,zlo)
# A = (xhi-xlo,0,0); B = (xy,yhi-ylo,0); C = (xz,yz,zhi-zlo)
# prism args = xlo xhi ylo yhi zlo zhi xy xz yz
region   box  prism 0 ${a11_box} 0 ${a22_box} 0 ${a33_box} ${a21_box} ${a31_box} ${a32_box} units box

create_box 1 box
create_atoms 1 box

# 把所有原子整体平移一个微量 ${shift}，使角原子 (0,0,0) 离开 xlo/ylo/zlo 面。
# 否则在 minimize / box-relax 时浮点噪声会把角原子推到盒子外侧（lambda<0），
# write_restart/read_restart 便会丢原子并报 "Did not assign all restart atoms correctly"。
# 这是周期性体相的纯刚性平移，能量/受力不变（对应 create.py 里的 INPLANE_SHIFT）。
displace_atoms all move ${shift} ${shift} ${shift} units box

# mass 已移至 general_mass.mod（经 general_potential.mod 统一 include）