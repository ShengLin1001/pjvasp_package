#!/home/yin/opt/bin/python3

import numpy as np 
from myvasp import vasp_func as vf 
import pandas as pd 
import shutil, os



def main():
    jobn, Etot, Eent, pres = vf.vasp_read_post_data()
    latoms2 = vf.get_list_of_outcar()

    try: 
        latoms2[0].get_magnetic_moments()
        plot_output(latoms2, jobn)
    except:
        print('==> no mag info.')




def plot_output(latoms2, jobn):
    try:
        shutil.rmtree('./y_post_mag/')
    except:
        print('==> no folder')
    os.mkdir('./y_post_mag/')


    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
     
    fig_wh = [3.15, 3]
    fig_subp = [1, 1]

    mag0_all = np.array([])
    
    for i in np.arange( len(jobn) ):
        magmom = latoms2[i].get_magnetic_moments()
        mag0 = latoms2[i].get_magnetic_moment() \
            / len(latoms2[i].get_positions())

        mag0_all = np.append(mag0_all, mag0)

        chem_sym = latoms2[i].get_chemical_symbols()
        elem_sym = pd.unique( chem_sym )
        
        atom_num = latoms2[i].get_atomic_numbers()
        elem_num = pd.unique( atom_num )
        nelem = len(elem_num)
        # print('elem_sym, elem_num, nelem:', \
        #     elem_sym, elem_num, nelem)

        ymax = np.ceil( max([1.5, abs(max(magmom, key=abs)) ]) )

        fig1, ax1 = vf.my_plot(fig_wh, fig_subp)

        for j in np.arange(nelem):
            mask = (atom_num == elem_num[j])
            yi = magmom[mask]
            temp = np.linspace(-0.2, 0.2, len(yi)+2 ) + j
            xi = temp[1:-1]
            ax1.plot( xi , yi, 'o' )

            ax1.text(j, ymax/4*5.0, '%.2f' %( yi.mean() ), \
                horizontalalignment='center' )
            ax1.text(j, ymax/4*4.5, '%.2f' %( yi.std() ), \
                horizontalalignment='center' )

            
        ax1.text(-0.5, ymax/4*5.0, 'mean:', \
            horizontalalignment='right' )
        ax1.text(-0.5, ymax/4*4.5, 'std:', \
            horizontalalignment='right' )
        ax1.text(-0.5, ymax/4*(-5.3), 'Supercell mean= %.2f' %(mag0), \
            horizontalalignment='left' )
            
        ax1.plot( [-0.5, nelem-0.5], [0, 0], '--k' )

        ax1.set_position([0.23, 0.15, 0.72, 0.70])

        ax1.set_xlim([-0.5, nelem-0.5])
        ax1.set_ylim([-ymax, ymax])
        ax1.set_xticks(np.arange(nelem))
        ax1.set_xticklabels(elem_sym)
        ax1.set_ylabel('Atomic magnetic moment ($\\mu_B$)')

        filename = './y_post_mag/y_post_mag_%s.pdf' %(jobn[i])
        plt.savefig(filename)
        plt.close('all')


    mag0_all = np.abs( mag0_all )
    f = open('./y_post_mag/y_post_mag_note.txt', 'w+')
    f.write('# VASP magnetic moment , supercell mean (mu_B/atom): \n' )
    
    f.write('\n%8s %8s %8s \n' \
        %('mean', 'max', 'min' ) )
    f.write('%8.2f %8.2f %8.2f \n\n' \
        %(mag0_all.mean(), mag0_all.max(), mag0_all.min() ) )
    
    for i in np.arange( len(jobn) ):
        f.write('%20s %8.2f \n' \
            %(jobn[i], mag0_all[i]  ) )
    f.close()



main()



