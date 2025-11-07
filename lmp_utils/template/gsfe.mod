############################ gsfe loop ############################


################# Check odd step #################
# Custom set parameters, for cubic plot
# v_step must be odd number

variable step equal 21
variable large_strain equal -2/3

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
variable strain  equal ${large_strain}*${factor}

# lammps cell
# xx 0  0
# xy yy 0
# xz yz zz
# Here we change yz from 0 to -2/3 * (yhi-ylo), to simulate the fcc(111) GSFE and hcp (0001) GSFE
variable deltayz equal  ${strain}*ly
change_box all yz delta ${deltayz}  units box

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
