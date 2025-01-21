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

    if jobn[0] != '00':
        sys.exit('==> ABORT! no reference state. ')
        
    latoms = vf.get_list_of_atoms()
    
    Asf = np.linalg.norm( \
        np.cross(latoms[0].cell[0, :], latoms[0].cell[1, :] ) )
    a11 = latoms[0].cell[0, 0]
    a22 = latoms[0].cell[1, 1]

    natoms = latoms[0].get_positions().shape[0]    
    E0bulk = Etot[0] / natoms


    dE, da3, dpos3 = check_constraints(Etot, latoms)
   
    gamma = dE/Asf *qe*1e23   
    sf, usf = find_sf_usf(gamma)

    from ase.formula import Formula
    str2 = latoms[0].get_chemical_formula()
    str2 = Formula(str2).format('latex')
    str_all = '%s\n$A =$%.4f $\\mathrm{\\AA}^2$' %(str2, Asf)  

    
    #=========================
    write_output(Asf, a11, a22, E0bulk, sf, usf, jobn, dE, gamma, da3)
    plot_GSFE(jobn, gamma, da3, dpos3, latoms, str_all)








def check_constraints(Etot, latoms):
    njobs = len(latoms)
    natoms = latoms[0].positions.shape[0]
    print('njobs, natoms:', njobs, natoms)

    dE = np.array([])
    da3 = np.zeros([1, 3])
    dpos3 = np.zeros([1, natoms])
    
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
        temp = np.linalg.norm(dlatt[0:2, :])
        if temp > 1e-10:
            print('dlatt:', dlatt)
            print('\n==> i, norm: {0}'.format([i, temp]) )
            sys.exit("==> ABORT: in-plane lattices relaxed. \n" )
    
        temp = dlatt[2,:].copy()
        da3 = np.vstack([ da3, temp[np.newaxis,:] ])
        

        # check pos
        dpos = latoms[i].positions - latoms[0].positions
        dposD = dpos @ np.linalg.inv(latoms[i].cell[:])
        dposD = dposD - np.around(dposD)  
        dpos = dposD @ latoms[i].cell[:]

        temp = np.linalg.norm(dpos[:,0:2])
        if temp > 1e-10:
            print('dpos', dpos)
            print('\n==> i, norm: {0}'.format([i, temp]) )
            sys.exit("==> ABORT: atoms show in-plane relaxation. \n" )
    
        temp = dpos[:, 2].copy()
        temp = temp - temp.mean()
        dpos3 = np.vstack([ dpos3, temp[np.newaxis, :] ])
        
    da3 = np.delete(da3, 0, 0)
    dpos3 = np.delete(dpos3, 0, 0)
    
    if (dE.shape[0] != njobs)  \
        or (da3.shape[0] != njobs) \
        or (da3.shape[1] != 3) \
        or (dpos3.shape[0] != njobs) \
        or (dpos3.shape[1] != natoms):
        sys.exit("==> ABORT: wrong dimensions. \n" )

    for i in np.arange(njobs):
        temp = np.abs(da3[i, 1]*da3[1, 0] - da3[1, 1]*da3[i, 0])
        if temp > 1e-10:
            print('\n==> i, norm: {0}'.format([i, temp]) )
            sys.exit("==> ABORT: slip is not along a line. \n" )
    
    return dE, da3, dpos3



def find_sf_usf(gamma):
    njobs = gamma.shape[0]
    print(njobs)
    sf = np.array([])
    usf = np.array([])
    for i in np.arange(1, njobs-1, 1):
        if (gamma[i] < gamma[i-1]) and (gamma[i] < gamma[i+1]):
            sf = np.append(sf, gamma[i])
        if (gamma[i] > gamma[i-1]) and (gamma[i] > gamma[i+1]):
            usf = np.append(usf, gamma[i])
    return sf, usf


   

def write_output(Asf, a11, a22, E0bulk, sf, usf, jobn, dE, gamma, da3):
    njobs = gamma.shape[0]
    print(njobs)
       
    f = open('y_post_gsfe.txt','w+')
    f.write('# VASP GSFE: \n' )
    f.write('# gamma = dE/Asf \n' )


    f.write('\n%16s %16s %16s %16s \n' \
        %('Asf (Ang^2)', 'a11 (Ang)', 'a22 (Ang)', 'E0_bulk (eV)' ) )
    f.write('%16.8f %16.8f %16.8f %16.8f \n\n' \
        %(Asf, a11, a22, E0bulk ))


    if sf.shape[0] > 0.5 :
        f.write('%20s: %16.8f \n' %('local min (mJ/m^2)', sf.min() ) )
    
    if usf.shape[0] > 0.5:
        f.write('%20s: %16.8f \n' %('local max (mJ/m^2)', usf.min() ) )
        
        
    f.write('\n%10s %10s %16s %10s %10s %10s %10s \n' \
        %('jobn', 'dE (eV)', 'gamma (mJ/m^2)', 
        'da31/a11', 'da32/a22', 'slip (Ang)',
        'da33 (Ang)'  
        ))
    
    for i in np.arange(njobs):
        f.write('%10s %10.4f %16.8f %10.4f %10.4f %10.4f %10.4f \n' \
            %(jobn[i], dE[i], gamma[i], 
            da3[i, 0]/a11,  da3[i, 1]/a22, np.linalg.norm(da3[i, 0:2]),
            da3[i, 2] 
            ))

    f.write(' \n')
    f.close() 





def plot_GSFE(jobn, gamma, da3, dpos3, latoms, str_all):
    njobs = gamma.shape[0]
    print(njobs)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig_wh = [3.15, 7]
    fig_subp = [3, 1]
    fig1, ax1 = vf.my_plot(fig_wh, fig_subp)

    xi = np.array([])
    disp = np.array([])
    for i in np.arange( njobs ):
        xi = np.append(xi, float(jobn[i]) )
        disp = np.append(disp, np.linalg.norm(da3[i, 0:2]))

    xi = xi / np.around( np.max([xi.max(), 10]) , -1)

    ax1[0].plot(xi, gamma, '-o')   
    ax1[1].plot(xi, da3[:,2], '-o')
      
    tau = np.diff(gamma) / np.diff(disp) *1e-2  #[GPa]
    x_tau = xi[0:-1].copy() + np.diff(xi)/2
    ax1[2].plot(x_tau, tau, '-o')
    ax1[2].plot([xi.min(), xi.max()], [0, 0], '--k')

    fig_pos  = np.array([0.23, 0.70, 0.70, 0.265])
    fig_dpos = np.array([0, -0.31, 0, 0])

    for i in np.arange(3):
        ax1[i].set_position(fig_pos + i*fig_dpos)
    
    ax1[-1].set_xlabel('Normalized slip vector')
    ax1[0].set_ylabel('GSFE (mJ/m$^2$)')
    ax1[1].set_ylabel('Inelastic normal displacement ($\\mathrm{\\AA}$)')
    ax1[2].set_ylabel('Shear stress $\\tau$ (GPa)')

    dxi = np.around( xi.max()/6, 1)
    ax1[-1].set_xticks( np.arange(0, xi.max()+dxi, dxi ) )


    ax1[2].text(0, tau.min()/2, \
        '$\\tau_\\mathrm{max} =$ %.1f GPa' %( tau.max() ))


    ax1[0].text( 0.3, gamma.max()*0.2, str_all )





    plt.savefig('y_post_gsfe.pdf')


    #=====================
    fig_wh = [5, 5]
    fig_subp = [1, 1]
    fig2, ax2 = vf.my_plot(fig_wh, fig_subp)

    xi = latoms[0].positions[:, 2]
    for i in np.arange(njobs):
        temp = np.array([ xi, dpos3[i, :] ])
        ind = np.argsort(temp[0, :])
        temp2 = temp[:,ind]

        ax2.plot(temp2[0, :], temp2[1, :], '-o', label=jobn[i])
    
    ax2.legend(loc='lower center', ncol=5, framealpha=0.4)
    ax2.set_xlabel('Atom positions in $x_3$ ($\\mathrm{\\AA}$)')
    ax2.set_ylabel('Displacement $u_3$ ($\\mathrm{\\AA}$)')
    ax2.set_position([0.17, 0.10, 0.78, 0.86])

    plt.savefig('y_post_gsfe.u3.pdf')



main()



