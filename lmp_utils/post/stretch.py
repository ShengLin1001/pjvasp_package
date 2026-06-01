import os
from mymetal.post.stretch import post_lammps_stretch


for lat in ['hcp', 'bcc', 'fcc']:
    dumpdir = os.path.join('./y_stretch', lat, 'dump')
    outputdir = os.path.join('./y_stretch', lat)
    #print(os.getcwd())
    post_lammps_stretch(post_data_file=os.path.join(dumpdir, 'y_post_data.txt'), 
                        refcontcar=os.path.join(dumpdir, 'data.full_relax'),
                        latoms_lammpstrj=os.path.join(dumpdir, 'movie.lammpstrj'),
                        save_fig_path=os.path.join(outputdir, 'p_post_stretch.pdf'),
                        save_txt_path=os.path.join(outputdir, 'p_post_stretch.txt')
                        )
