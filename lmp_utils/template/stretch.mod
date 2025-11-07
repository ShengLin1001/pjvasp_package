############################ stretch loop ############################

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

# To stretch in-plane and/or out-of-plane lattice
variable deltaxx equal ${strain}*lx
variable deltayy equal ${strain}*ly
variable deltazz equal ${strain}*lz

variable deltaxy equal ${strain}*xy
variable deltaxz equal ${strain}*xz
variable deltayz equal ${strain}*yz

# iso expand/compress box
change_box all x delta 0 ${deltaxx} y delta 0 ${deltayy} z delta 0 ${deltazz} &
               xy delta ${deltaxy} xz delta ${deltaxz} yz delta ${deltayz} remap units box

# For hcp, relax c axis only
# For fcc/bcc, relax none of the box lengths
include stretch_constrained_relax.mod

include general_structural_info.mod

variable jobn string "$(v_stretch:%16.12e)"
include general_output.mod

next i
jump SELF loop
