############################ stretch loop ############################


################# Check odd step #################
# Custom set parameters, for cubic plot
# v_step must be odd number

variable step equal 101
variable large_strain equal 0.003

variable check_even equal ${step}%2

if "${check_even} < 1e-10" then &
    "variable temp equal ${step}+1" &
else &
    "variable temp equal ${step}" &

variable step equal ${temp}

################## Check end ##################


################### Stretch loop ###################
label loop
clear
read_restart ../dump/restart.equil
include general_potential.mod

# Get stretch factor v_stretch
variable i  loop ${step}                  # 1, 2, ..., step
variable id equal ${i}-(${step}-1)/2-1    # -50, -49, ..., 0, 1, ..., 50.
variable factor  equal ${id}/((${step}-1)/2)  # -1, -49/50, ..., 0, 1/50, ..., 1.
variable strain  equal ${large_strain}*${factor}
variable stretch equal ${strain}+1


################## case ##################
# dirn = 
# 1: y_cij_energy_c11
# 2: y_cij_energy_c12
# 3: y_cij_energy_c13
# 4: y_cij_energy_c33
# 5: y_cij_energy_c44

variable deltaxx equal ${strain}*lx
variable deltayy equal ${strain}*ly
variable deltazz equal ${strain}*lz
variable deltaxy equal ${strain}*xy
variable deltaxz equal ${strain}*xz
variable deltayz equal ${strain}*yz
# lammps cell
# xx 0  0
# xy yy 0
# xz yz zz
if "${dirn} == 1" then &
    "change_box all x delta 0 ${deltaxx} xy delta ${deltaxy} xz delta ${deltaxz} remap units box"    &
    "variable filename string y_cij_energy_c11" &
elif "${dirn} == 2"  &
    "change_box all x delta 0 ${deltaxx} xy delta ${deltaxy} xz delta ${deltaxz} &
                    y delta 0 ${deltayy} yz delta ${deltayz}                     remap units box"    &
    "variable filename string y_cij_energy_c12" &
elif "${dirn} == 3"  &
    "change_box all x delta 0 ${deltaxx} xy delta ${deltaxy} xz delta ${deltaxz} &
                    z delta 0 ${deltazz}                                         remap units box"    &
    "variable filename string y_cij_energy_c13" &
elif "${dirn} == 4"  &
    "change_box all z delta 0 ${deltazz}                                         remap units box"    &
    "variable filename string y_cij_energy_c33" &
elif "${dirn} == 5"  &
    "variable deltayz equal ${strain}*lz"      &
    "change_box all yz delta ${deltayz}                                          remap units box"    &
    "variable filename string y_cij_energy_c44" &
else                 &
    "print 'Error: dirn wrong!'"


################## case end ##################

min_style   cg
minimize ${etol} ${ftol} ${maxiter} ${maxeval}

variable jobn  string "$(v_stretch:%16.12e)"
include Cij_energy_output.mod

next i
jump SELF loop
