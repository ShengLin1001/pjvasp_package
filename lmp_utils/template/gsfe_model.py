from mymetal.build.bulk.create import create_fcc_111, create_hcp_basal, create_hcp_prism1
from ase.io import write

a = aa_template
c = cc_template

latnum = lat_template

if latnum == 1:  # hcp
    atoms = create_hcp_basal(a=a, c=c, size = (1, 1, 12))
elif latnum == 2:  # fcc
    atoms = create_fcc_111(a=a, size = (1, 1, 8))

write('data.ini', atoms, format='lammps-data')
