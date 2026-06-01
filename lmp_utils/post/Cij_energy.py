from mymetal.post.Cij_energy import post_lammps_Cij_energy
import os
    
for lat in ['hcp', 'bcc', 'fcc']:
    workdir = os.path.join('./y_Cij_energy', lat)
    dumpdir = os.path.join('./y_Cij_energy', lat, 'dump') # Contains deformation subdirs

    post_lammps_Cij_energy(dir=dumpdir,
                        refcontcar=os.path.join(dumpdir, 'data.full_relax'),
                        save_fig_path=os.path.join(workdir, 'y_post_cij_energy.pdf'),
                        save_txt_path=os.path.join(workdir, 'y_post_cij_energy.txt'),)
