# Lattice parameters
# initial lattice constant
variable aa equal aa_template
variable lat equal lat_template   # 1 - "hcp", 2 - "fcc", 3 - "bcc"

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
mass 1 196.97