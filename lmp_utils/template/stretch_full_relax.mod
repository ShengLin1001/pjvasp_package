
if "${lat} == 1" then &
    "jump SELF hcp"     &
elif "${lat} == 2"    &
    "jump SELF bcc_fcc"     &
elif "${lat} == 3"    &
    "jump SELF bcc_fcc"     &
else &
    "print 'Unknown lattice type: ${lat}'"

label bcc_fcc

# scaleyz, scalexz, scalexy = yes, to keep the angles fixed
fix br all box/relax  iso 0 fixedpoint 0 0 0

jump SELF end_if


label hcp

# scaleyz, scalexz, scalexy = yes, to keep the angles fixed
fix br all box/relax  x 0 y 0 z 0 couple xy fixedpoint 0 0 0

jump SELF end_if

label end_if

min_style   cg
minimize ${etol} ${ftol} ${maxiter} ${maxeval}
unfix       br
write_data    ../dump/data.full_relax
write_restart ../dump/restart.equil