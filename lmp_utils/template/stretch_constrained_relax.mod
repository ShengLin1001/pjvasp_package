
if "${lat} == 1" then &
    "jump SELF hcp"     &
elif "${lat} == 2"    &
    "jump SELF bcc_fcc"     &
elif "${lat} == 3"    &
    "jump SELF bcc_fcc"     &
else &
    "print 'Unknown lattice type: ${lat}'"

label bcc_fcc

#fix br all box/relax  iso 0 fixedpoint 0 0 0

jump SELF end_if


label hcp

# scaleyz, scalexz, scalexy = yes, to keep the angles fixed
fix br all box/relax  z 0 fixedpoint 0 0 0

jump SELF end_if

label end_if

min_style   cg
minimize ${etol} ${ftol} ${maxiter} ${maxeval}

if "${lat} == 1" then "unfix br"