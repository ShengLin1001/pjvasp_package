import os
from mymetal.post.gsfe import post_gsfe

for lat in ['fcc', 'hcp']:
    
    if lat == 'fcc':
        lat_types = ['FCC_111', 'FCC_100']
    elif lat == 'hcp':
        lat_types = ['HCP_basal', 'HCP_prism1w', 'HCP_pyr1w', 'HCP_pyr2']

    for gsfe_type in lat_types:
        workdir = os.path.join('./y_gsfe', lat, gsfe_type)
        dumpdir = os.path.join('./y_gsfe', lat, gsfe_type, 'dump')

        post_gsfe(
                  post_data_path=os.path.join(dumpdir, 'y_post_data.txt'),
                    latoms_path=os.path.join(dumpdir, 'movie.lammpstrj'),
                    latoms_format='lammps-dump-text',
                    save_fig_path=os.path.join(workdir, 'y_post_gsfe.pdf'),
                    save_u3_fig_path=os.path.join(workdir, 'y_post_gsfe.u3.pdf'),
                    save_txt_path=os.path.join(workdir, 'y_post_gsfe.txt'),)

