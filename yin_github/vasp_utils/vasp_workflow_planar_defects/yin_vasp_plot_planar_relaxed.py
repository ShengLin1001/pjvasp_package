#!/public1/home/a8s000114/opt/bin/python3


import numpy as np
from myvasp import vasp_func as vf 
import sys


def main():
    qe = vf.phy_const('qe')
    
    jobn, Etot, Eent, pres = vf.vasp_read_post_data()
    njobs = len(jobn)  # number of jobs

    if njobs < 1.5:
        sys.exit('==> ABORT! more structures needed. ')

    ibulk=-1
    for i in np.arange(njobs):
        if jobn[i] == 'bulk':
            ibulk = i 
            break

    if ibulk == -1:
        sys.exit('==> ABORT! no reference bulk state. ')
        
    latoms = vf.get_list_of_atoms()
    
    Asf = np.linalg.norm( \
        np.cross(latoms[ibulk].cell[0, :], latoms[ibulk].cell[1, :] ) )
    a11 = latoms[ibulk].cell[0, 0]
    a22 = latoms[ibulk].cell[1, 1]
    
    natoms = latoms[ibulk].get_positions().shape[0]
    E0bulk = Etot[ibulk]/natoms
    V0bulk = latoms[ibulk].get_volume()/natoms

    if np.abs( Asf-a11*a22 ) > 1e-10:
        sys.exit('ABORT: wrong Asf. ')


    dE, da3, dpos_all = check_constraints(Etot, latoms, ibulk)
   
    gamma = dE/Asf *qe*1e23   #[mJ/m^2]
    
    #=========================
    write_output(Asf, a11, a22, E0bulk, V0bulk, jobn, dE, gamma, da3)
    plot_output(jobn, latoms, dpos_all, gamma, Asf, ibulk)







def check_constraints(Etot, latoms, ibulk):
    njobs = len(latoms)
    natoms = latoms[ibulk].positions.shape[0]
    print('njobs, natoms:', njobs, natoms)

    dE = np.array([])
    da3 = np.zeros([1, 3])
    dpos_all = np.zeros([1, natoms, 3])
    
    for i in np.arange(njobs):
        dE = np.append(dE, Etot[i]-Etot[ibulk])
    

        # check elem
        if latoms[i].get_chemical_formula() \
            != latoms[ibulk].get_chemical_formula():
            sys.exit('ABORT: wrong chemical formula. ')

        temp = latoms[i].get_atomic_numbers() \
            - latoms[ibulk].get_atomic_numbers()
        vf.confirm_0( temp )
        

        # check latt 
        dlatt = latoms[i].cell[:] - latoms[ibulk].cell[:]

        temp = dlatt[0:2, :].copy()
        temp = np.linalg.norm(temp)
        if temp > 1e-10:
            print('dlatt:', dlatt)
            print('\n==> i, norm: {0}'.format([i, temp]) )
            sys.exit("==> ABORT: in-plane lattices changed. \n" )
    

        temp = dlatt[2, :].copy()
        da3 = np.vstack([ da3, temp[np.newaxis, :] ])
        

        # check pos
        dpos = latoms[i].positions - latoms[ibulk].positions
        dposD = dpos @ np.linalg.inv(latoms[i].cell[:])
        dposD = dposD - np.around(dposD)  
        dpos = dposD @ latoms[i].cell[:]

        temp = dpos.copy()
        for j in np.arange(3):
            temp[:,j] = temp[:,j] - temp[:,j].mean()
        dpos_all = np.vstack([ dpos_all, temp[np.newaxis, :] ])
        
          
    da3 = np.delete(da3, 0, 0)
    dpos_all = np.delete(dpos_all, 0, 0)
    
    if (dE.shape[0] != njobs)  \
        or (da3.shape[0] != njobs)  \
        or (dpos_all.shape[0] != njobs) \
        or (dpos_all.shape[1] != natoms) \
        or (dpos_all.shape[2] != 3) :
        sys.exit("==> ABORT: wrong dimensions. \n" )
   
    return dE, da3, dpos_all



   

def write_output(Asf, a11, a22, E0bulk, V0bulk, jobn, dE, gamma, da3):
    njobs = gamma.shape[0]
    print('njobs:', njobs)
       
    f = open('y_post_planar_relaxed.txt','w+')
    f.write('# VASP relaxed planar defects: \n' )
    f.write('# gamma = dE/Asf \n' )


    f.write('\n%16s %16s %16s \n' \
        %('Asf (Ang^2)', 'a11 (Ang)', 'a22 (Ang)' ) )
    f.write('%16.8f %16.8f %16.8f \n' \
        %(Asf, a11, a22 ))

    f.write('\n%16s %16s %16s %16s \n' \
        %('E0_bulk (eV)', 'V0_bulk (Ang^3)', 'a0_fcc', 'a0_bcc' ) )
    f.write('%16.8f %16.8f %16.8f %16.8f \n' \
        %(E0bulk, V0bulk, (V0bulk*4)**(1/3), (V0bulk*2)**(1/3) ))


    f.write('\n%10s %10s %16s %10s %10s %10s %10s %10s \n' \
        %('jobname', 'dE (eV)', 'gamma (mJ/m^2)', 'gamma/2', 
          'da31/a11', 'da32/a22', 'slip (Ang)', 
          'da33 (Ang)'
          ))
    
    for i in np.arange(njobs):
        f.write('%10s %10.4f %16.8f %10.4f %10.4f %10.4f %10.4f %10.4f \n' \
            %(jobn[i], dE[i], gamma[i], gamma[i]/2,
              da3[i, 0]/a11, da3[i, 1]/a22, np.linalg.norm(da3[i, 0:2]),  
              da3[i, 2] 
              ))

    f.write(' \n')
    f.close() 





def plot_output(jobn, latoms, dpos_all, gamma, Asf, ibulk):
    from ase.formula import Formula
    
    njobs = len(jobn)
    print('njobs:', njobs)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig_wh = [3.15, 3]
    fig_subp = [1, 1]

    xi = latoms[ibulk].positions[:, 2].copy()
    
    for i in np.arange(njobs):
        
        if np.linalg.norm( dpos_all[i, :] ) > 1e-10:
            fig1, ax1 = vf.my_plot(fig_wh, fig_subp)
            
            temp = np.hstack([ dpos_all[i, :], xi[:, np.newaxis] ])
            ind = np.argsort(temp[:, -1])
            temp = temp[ind, :]
    
            ax1.plot(temp[:, -1], temp[:, 0], '-s', label='$u_1$'  )
            ax1.plot(temp[:, -1], temp[:, 1], '-o', label='$u_2$'  )
            ax1.plot(temp[:, -1], temp[:, 2], '-^', label='$u_3$'  )

            ax1.legend(loc='lower center', ncol=3, framealpha=0.4)
    
            ax1.set_xlabel('Atom positions in $x_3$ ($\\mathrm{\\AA}$)')
            ax1.set_ylabel('Displacements $u_i$ ($\\mathrm{\\AA}$)')
            ax1.set_position([0.25, 0.16, 0.7, 0.76])


            if jobn[i] == 'ssf':
                str1 = '$\\gamma_\\mathrm{ssf} =$ %.0f mJ/m$^2$'  %(gamma[i])
            
            elif jobn[i] == 'usf':
                str1 = '$\\gamma_\\mathrm{usf} =$ %.0f mJ/m$^2$' %(gamma[i])
            
            elif jobn[i] == 'surf':
                str1 = '$\\gamma_\\mathrm{surf} =$ %.0f mJ/m$^2$' %(gamma[i]/2)

            elif jobn[i] == 'twin':
                str1 = '$\\gamma_\\mathrm{twin} =$ %.0f mJ/m$^2$' %(gamma[i]/2)
            
            else:
                str1 = '$\\Delta E / A =$ %.0f mJ/m$^2$'          %(gamma[i])
            
            str2 = latoms[ibulk].get_chemical_formula()
            str2 = Formula(str2).format('latex')
            
            str_all = '%s\n$A =$%.4f $\\mathrm{\\AA}^2$\n%s' \
                %(str2, Asf, str1)  
        
            ax1.text( xi.max()*0.2, dpos_all[i, :].max()*0.6, str_all )

            filename = 'y_post_planar_relaxed.%s.pdf' %(jobn[i])

            plt.savefig(filename)



main()



