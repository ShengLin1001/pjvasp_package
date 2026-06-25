from mymetal.build.bulk.gsfe import create_gsfe_model
from ase.io import write

a = aa_template
c = cc_template

latnum = lat_template

gsfe_type = "gsfe_type_template"

atoms = create_gsfe_model(gsfe_type=gsfe_type, a=a, c=c)

write('data.ini', atoms, format='lammps-data')
