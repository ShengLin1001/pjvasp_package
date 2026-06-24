# Lattice parameters
# initial lattice constant
variable aa equal aa_template
variable lat equal lat_template   # 1 - "hcp", 2 - "fcc", 3 - "bcc"
variable shift equal 1e-5

print "Lattice type: ${lat}"

if "${lat} == 1" then &
    "jump SELF hcp"     &
elif "${lat} == 2"    &
    "jump SELF fcc"     &
elif "${lat} == 3"    &
    "jump SELF bcc"     &
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