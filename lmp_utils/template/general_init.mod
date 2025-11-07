# Uncomment one of these blocks, depending on what units
# you are using in LAMMPS and for output

# metal units, elastic constants in GPa
units		metal
variable cfac equal 1.0e-4
variable cunits string GPa

# Define minimization parameters
variable etol equal 0.0 
variable ftol equal 1.0e-8
variable maxiter equal 30000
variable maxeval equal 10000

# generate settings
dimension       3
atom_style      atomic
#timestep    1.0
boundary	p p p
