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

    if jobn[0] != '0.00':
        sys.exit('==> ABORT! no reference state. ')
        
    latoms = vf.get_list_of_atoms()
    
    Asf = np.linalg.norm( \
        np.cross(latoms[0].cell[0, :], latoms[0].cell[1, :] ) )

    a11 = latoms[0].cell[0, 0]
    a22 = latoms[0].cell[1, 1]
    
    natoms = latoms[0].get_positions().shape[0]    
    E0bulk = Etot[0] / natoms


    dE, da33 = check_constraints(Etot, latoms)
   
    # check jobname
    k = np.array([])
    for i in np.arange(len(jobn)):
        k = np.append(k, float(jobn[i]) )
     
    if np.linalg.norm( k - da33 ) > 1e-10:
        sys.exit('==> ABORT. wrong jobname. ')


    gamma = dE/Asf *qe*1e23   #[mJ/m^2]
   
    from ase.formula import Formula
    str2 = latoms[0].get_chemical_formula()
    str2 = Formula(str2).format('latex')
    str_all = '%s\n$A =$%.4f $\\mathrm{\\AA}^2$' %(str2, Asf)

 
    #=========================
    write_output(Asf, a11, a22, E0bulk, jobn, dE, gamma, da33)
    plot_output(gamma, da33, str_all)







def check_constraints(Etot, latoms):
    njobs = len(latoms)
    natoms = latoms[0].positions.shape[0]
    print('njobs, natoms:', njobs, natoms)

    dE = np.array([])
    da33 = np.array([])

    for i in np.arange(njobs):
        dE = np.append(dE, Etot[i]-Etot[0])
    

    # check elem
        if latoms[i].get_chemical_formula() \
            != latoms[0].get_chemical_formula():
            sys.exit('ABORT: wrong chemical formula. ')

        temp = latoms[i].get_atomic_numbers() \
            - latoms[0].get_atomic_numbers()
        vf.confirm_0( temp )


        # check latt 
        dlatt = latoms[i].cell[:] - latoms[0].cell[:]

        temp  = dlatt[0:2, :].copy()
        temp2 = dlatt[2, 0:2].copy()

        temp  = np.linalg.norm(temp)
        temp2 = np.linalg.norm(temp2)

        if temp > 1e-10 or temp2 > 1e-10:
            print('dlatt:', dlatt)
            print(i, temp, temp2)
            sys.exit("==> ABORT: lattices wrong. \n" )
    
        da33 = np.append( da33, dlatt[2,2] )

        
        # check pos
        dpos = latoms[i].positions - latoms[0].positions
        dposD = dpos @ np.linalg.inv(latoms[i].cell[:])
        dposD = dposD - np.around(dposD)  
        dpos = dposD @ latoms[i].cell[:]

        temp = np.linalg.norm(dpos)
        if temp > 1e-10:
            print('dpos:', dpos)
            print('\n==> i, norm: {0}'.format([i, temp]) )
            sys.exit("==> ABORT: positions changed. \n" )
    
    
    if (dE.shape[0] != njobs):
        sys.exit("==> ABORT: wrong dimensions. \n" )
   
    return dE, da33




   

def write_output(Asf, a11, a22, E0bulk, jobn, dE, gamma, da33):
    njobs = gamma.shape[0]
    print(njobs)
       
    f = open('y_post_decohesion.txt','w+')
    f.write('# VASP decohesion: \n' )
    f.write('# gamma = dE/Asf \n' )

    f.write('\n%16s %16s %16s %16s \n' \
        %('Asf (Ang^2)', 'a11 (Ang)', 'a22 (Ang)', 'E0_bulk (eV)' ) )
    f.write('%16.8f %16.8f %16.8f %16.8f \n' \
        %(Asf, a11, a22, E0bulk ))

        
    f.write('\n%10s %10s %16s %10s %10s \n' \
        %('jobname', 'dE (eV)', 'gamma (mJ/m^2)', 'gamma/2', 
          'da33 (Ang)'  ))
    
    for i in np.arange(njobs):
        f.write('%10s %10.4f %16.8f %10.4f %10.4f \n' \
            %(jobn[i], dE[i], gamma[i], gamma[i]/2,
            da33[i] ))

    f.write(' \n')
    f.close() 





def plot_output(gamma, da33, str_all):
    njobs = len(gamma)
    print(njobs)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig_wh = [3.15, 5]
    fig_subp = [2, 1]
    fig1, ax1 = vf.my_plot(fig_wh, fig_subp)

    disp = da33.copy()
    xi = da33.copy()

    ax1[0].plot(xi, gamma, '-o')   
      
    tau = np.diff(gamma) / np.diff(disp) *1e-2  #[GPa]
    x_tau = xi[0:-1].copy() + np.diff(xi)/2
    ax1[1].plot(x_tau, tau, '-o')


    fig_pos  = np.array([0.23, 0.57, 0.70, 0.40])
    fig_dpos = np.array([0, -0.46, 0, 0])

    ax1[0].set_position(fig_pos)
    ax1[1].set_position(fig_pos + fig_dpos  )

    ax1[-1].set_xlabel('Vacuum layer thickness ($\\mathrm{\\AA}$)')
    ax1[0].set_ylabel('Decohesion energy (mJ/m$^2$)')
    ax1[1].set_ylabel('Tensile stress $\\sigma$ (GPa)')

    ax1[0].text(4, gamma.max()/2, \
        '$\\gamma_\\mathrm{surf}^\\mathrm{unrelaxed}$ \n= %.0f mJ/m$^2$' %( gamma[-1]/2 ))

    ax1[1].text(4, tau.max()/2, \
        '$\\sigma_\\mathrm{max}$ \n= %.1f GPa' %( tau.max() ))


    ax1[0].text( 2, gamma.max()*0.2, str_all )


    plt.savefig('y_post_decohesion.pdf')


       
    


main()



