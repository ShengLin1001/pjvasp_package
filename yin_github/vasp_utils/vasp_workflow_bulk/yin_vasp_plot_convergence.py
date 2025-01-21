
#!/home/yin/opt/bin/python3


import numpy as np
from myvasp import vasp_func as vf 
import os


def main():
   
    atoms_ref = vf.my_read_vasp('../y_full_relax/CONTCAR')  # for strain 
    natoms = atoms_ref.get_positions().shape[0]


    dirlist=[
        'y_convergence_encut', 
        'y_convergence_kp', 
        'y_convergence_sigma', 
    ]

    ljobn = []   # list of data 
    lEtot = [] 
    lEent = [] 


    for i in np.arange(len(dirlist)):
        dirn = dirlist[i]
        print(dirn) 

        os.chdir(dirn)
        jobn, Etot, Eent, pres = vf.vasp_read_post_data()
        os.chdir('..')

        x  = np.array([])     
        for j in np.arange(len(jobn)):
            temp = float(jobn[j])
            x = np.append(x, temp )

        ljobn.append(x) 
        lEtot.append(Etot/natoms) 
        lEent.append(Eent/natoms)  


    plot_convergence(ljobn, lEtot, lEent) 








def plot_convergence(ljobn, lEtot, lEent):

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt


    fig_wh = [7, 7]
    fig_subp = [3, 2]
    fig1, ax1 = vf.my_plot(fig_wh, fig_subp, fig_sharex=False)

    fig_pos  = np.array([0.12, 0.73, 0.35, 0.25])
    fig_dpos = np.array([0, -1/3, 0, 0])
    fig_dpos2 = np.array([0.5, 0, 0, 0])
 
    for i in np.arange( len(ljobn) ):
        ax1[i, 0].set_position(fig_pos + i* fig_dpos)
        ax1[i, 1].set_position(fig_pos + i* fig_dpos + fig_dpos2)



    for i in np.arange(2):
        ax1[i, 0].plot( ljobn[i], lEtot[i], '-o')
        ax1[i, 1].plot( ljobn[i][0:-1], np.diff( lEtot[i] )*1e3, '-o')


    ax1[2, 0].plot( ljobn[2], lEent[2], '-s')
    ax1[2, 1].plot( ljobn[2][0:-1], np.diff( lEent[2] )*1e3, '-s')


    for j in np.arange(2):
        ax1[0, j].set_xlabel('ENCUT')
        ax1[1, j].set_xlabel('KP')
        ax1[2, j].set_xlabel('SIGMA')

        ax1[0, j].set_xlim([0, 1000])
        ax1[1, j].set_xlim([0, 150])
        ax1[2, j].set_xlim([0, 0.5])

    for i in np.arange(3):
        ax1[i, 1].set_ylim([-1, 1])


    ax1[0, 0].set_ylabel('Energy (eV/atom)')
    ax1[1, 0].set_ylabel('Energy (eV/atom)')
    ax1[2, 0].set_ylabel('EENTRO (eV/atom)')

    ax1[0, 1].set_ylabel('$\\Delta_{(i+1)-(i)}$ (meV/atom)')
    ax1[1, 1].set_ylabel('$\\Delta_{(i+1)-(i)}$ (meV/atom)')
    ax1[2, 1].set_ylabel('$\\Delta_{(i+1)-(i)}$ (meV/atom)')


    plt.savefig('y_post_convergence.pdf')
    plt.close('all')




main() 




