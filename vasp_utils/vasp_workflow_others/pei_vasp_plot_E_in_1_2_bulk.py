import os
from ase.io import write, read
from myvasp import vasp_func as vf 
from mymetal.post.general import my_read_y_dir_contcar
from mymetal.post.E_in_1_2_bulk import post_E_in_1_2_bulk
# Pre-pare data

################### This part maybe need to be changed

# e. g. "****/20231215_strain_energy_fcc_hcp/more_strain_uniaxial_biaxial_bulk_fcc_hcp/fcc"
# Already cd the workflow root
os.chdir("..")
myroot = os.getcwd()
workflow_root = os.path.join(myroot, "y_E_in_1_2_bulk")
post_data_path = os.path.join(workflow_root, "y_post_data.txt")
y_dir_path = os.path.join(workflow_root, "y_dir")
refcontcar_path = os.path.join(myroot, "y_full_relax/CONTCAR")
latoms_path = os.path.join(workflow_root, "movie.xyz")
save_fig_path = os.path.join(workflow_root, "p_post_E_in_1_2_bulk.pdf")
save_fig_path2 = os.path.join(workflow_root, "p_post_E_in_1_2_bulk2.pdf")
save_txt_path = os.path.join(workflow_root, "p_post_E_in_1_2_bulk.txt")
################### To here

os.chdir(workflow_root)
os.system("pei_vasp_univ_post")

jobn, Etot, Eent, pres = vf.vasp_read_post_data(post_data_path)
latoms = my_read_y_dir_contcar(dir = y_dir_path, post_data_file = post_data_path)
atoms_ref = read(refcontcar_path)
# For viewing
write(latoms_path, latoms, format='extxyz')

post_E_in_1_2_bulk(jobn = jobn, Etot = Etot,
                   atoms_ref = atoms_ref,
                   latoms = latoms,
                   save_fig_path = save_fig_path, 
                   save_fig_path2 = save_fig_path2,
                   save_txt_path = save_txt_path
                   )