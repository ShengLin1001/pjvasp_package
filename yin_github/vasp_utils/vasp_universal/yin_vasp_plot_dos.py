
import numpy as np 
from myvasp import vasp_func as vf 
import shutil, os 


def main():
    try:
        shutil.rmtree('./y_post_dos/')
    except:
        print('==> no folder')
    os.mkdir('./y_post_dos/')

    jobn, Etot, Eent, pres = vf.vasp_read_post_data()
    njobs = len(jobn)

    for i in np.arange(njobs):
        fname = 'y_dir/%s/DOSCAR' %(jobn[i])

        if os.path.exists(fname):
            atoms_dos = vf.my_read_doscar(fname=fname)
            atoms_dos.plot_dos()
            
            file_new = './y_post_dos/y_post_dos_%s.pdf'  %( jobn[i] ) 
            shutil.move('y_post_dos.pdf', file_new)



main()



