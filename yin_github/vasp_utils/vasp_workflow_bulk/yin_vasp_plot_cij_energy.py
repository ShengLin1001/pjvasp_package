#!/home/yin/opt/bin/python3


import numpy as np
from myvasp import vasp_func as vf 
import os


def main():
  
    try:
        atoms_ref = vf.my_read_vasp('../y_full_relax/CONTCAR')  # for strain 
        print('reading ref CONTCAR...') 
    except:
        atoms_ref = vf.my_read_vasp('../y_full_relax/POSCAR')  # for strain
        print('reading ref POSCAR...')

 
    dirlist=[
        'y_cij_energy_c11', 
        'y_cij_energy_c12', 
        'y_cij_energy_c13', 
        'y_cij_energy_c33', 
        'y_cij_energy_c44', 
    ]

    ldata = []   # list of data 

    for i in np.arange(len(dirlist)):
        dirn = dirlist[i]
        print(dirn) 
        data = read_deform_data(dirn, atoms_ref)
        ldata.append(data) 

    check_ldata(ldata)

    plot_cij_energy(ldata)
   







def read_deform_data(dirn, atoms_ref):

    x  = np.array([])      # applied strain, with respect to the ref state 
    e_ss_V = np.zeros(6)   # strain components 
    
    os.chdir(dirn)
    jobn, Etot, Eent, pres = vf.vasp_read_post_data()
    latoms = vf.get_list_of_atoms()

    iref = -1    # id of the ref state 
    
    for i in np.arange(len(jobn)):
        temp = float(jobn[i])-1
        x = np.append(x, temp )
        
        if np.abs(temp)<1e-10:
            iref = i 

        temp = vf.calc_strain(atoms_ref.get_cell()[:], latoms[i].get_cell()[:])
        e_ss_V = np.vstack([ e_ss_V, temp])
    

    if iref < 0:
        os.exit('ABORT: no ref state. ')    
    
    e_ss_V = np.delete(e_ss_V, 0, 0)
    np.savetxt('y_post_calc_strain_all.txt', e_ss_V )
    os.chdir('..')


    # energy density
    Etot = Etot - Etot[iref]
    ed   = Etot / latoms[iref].get_volume()  * vf.phy_const('qe')*1e21    
        
    # stress     
    s = pres *(-0.1)   # GPa

    temp = s.copy()
    s[:,3] = temp[:,4]
    s[:,4] = temp[:,5]
    s[:,5] = temp[:,3]
    
      
    data = {}
    data['e']  = x
    data['ed'] = ed
    data['s']  = s
    data['iref']  = iref 
    return data 

    




def check_ldata(ldata):

    iref0 = ldata[0]['iref']

    for i in np.arange(len(ldata)):
        iref1 = ldata[i]['iref']
        
        vf.confirm_0( ldata[0]['e'][iref0]  - ldata[i]['e'][iref1]  )
        vf.confirm_0( ldata[0]['ed'][iref0] - ldata[i]['ed'][iref1] )
        vf.confirm_0( ldata[0]['s'][iref0]  - ldata[i]['s'][iref1]  )










def plot_cij_energy(ldata):

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    xi = np.linspace(-0.003, 0.003, 1000)

    fig_wh = [7, 10]
    fig_subp = [5, 2]
    fig1, ax1 = vf.my_plot(fig_wh, fig_subp)

    fig_pos  = np.array([0.13, 0.8, 0.35, 0.16])
    fig_dpos = np.array([0, -0.19, 0, 0])
    fig_dpos2 = np.array([0.5, 0, 0, 0])
 
    for i in np.arange(5):
        ax1[i, 0].set_position(fig_pos + i* fig_dpos)
        ax1[i, 1].set_position(fig_pos + i* fig_dpos + fig_dpos2)



    # plot energy-strain 

    p0_tot = np.array([])     # for Cij
    e0_tot = np.array([])     # to check shift 

    for i in np.arange(len(ldata)):
        e  = ldata[i]['e']
        ed = ldata[i]['ed'] 
        s  = ldata[i]['s'] 

        ax1[i, 0].plot(e, ed, 'o')

        param = np.polyfit(e, ed, 2)
        fun = np.poly1d(param)
        yi = fun(xi)
        p0_tot = np.append(p0_tot, param[0]) 

        e0 = -param[1]/(2*param[0])
        e0_tot = np.append(e0_tot, e0 ) 

        ax1[i, 0].plot(xi, yi, '-C1')

        str1 = '$\\epsilon_0$ = %.6f' %(e0)
        ymin, ymax = ax1[i, 0].get_ylim()
        ax1[i, 0].text(0, ymax*0.85, str1)

    C11 = p0_tot[0]*2
    C12 = p0_tot[1]-C11
    C33 = p0_tot[3]*2
    C13 = p0_tot[2]-(C11+C33)/2
    C44 = p0_tot[4]*2



    # plot stress-strain 

    lcolor=['C0','C1','C2','C3','C4','C5']

    llegend=['$\\sigma_{xx}$','$\\sigma_{yy}$','$\\sigma_{zz}$','$\\sigma_{yz}$','$\\sigma_{xz}$','$\\sigma_{xy}$']

    lslope = np.array([
        [ C11    , C12    , C13    , 0, 0, 0 ],          # c11
        [ C11+C12, C12+C11, C13*2  , 0, 0, 0 ],          # c12
        [ C11+C13, C12+C13, C13+C33, 0, 0, 0 ],          # c13
        [ C13    , C13    , C33    , 0, 0, 0 ],          # c33
        [ 0      , 0      , 0      , C44,  0,  0 ],      # c44
    ])
 
    for i in np.arange(len(ldata)):
        e    = ldata[i]['e']
        ed   = ldata[i]['ed'] 
        s    = ldata[i]['s'] 
        iref = ldata[i]['iref'] 

        for j in np.arange(6):
            ax1[i, 1].plot(e, s[:, j], 'o', color=lcolor[j], label = llegend[j] )
            ax1[i, 1].plot(xi, xi * lslope[i,j] + s[iref,j], '-', color=lcolor[j] )

    ax1[0,1].legend(loc='upper left', ncol=2, framealpha=0.4)



    for i in np.arange(len(ldata)):
        for j in np.arange(2):
            ax1[i,j].set_xticks(np.arange(-0.002, 0.004, 0.002))
            ax1[4,j].set_xlabel('Strain')

        ax1[i,0].set_ylabel('Elastic energy density (GPa)')
        ax1[i,1].set_ylabel('DFT stress (GPa)')




    str1 = '$C_{11}$, $C_{12}$, $C_{13}$, $C_{33}$, $C_{44}$ (GPa): \n%.2f,   %.2f,   %.2f,   %.2f,   %.2f ' \
        %( C11, C12, C13, C33, C44 )
    ymin, ymax = ax1[0,0].get_ylim()
    ax1[0,0].text(-0.003, ymax*1.05, str1 )



    plt.savefig('y_post_cij_energy.pdf')
    plt.close('all')




    write_cij_energy( np.array([ C11, C12, C13, C33, C44]) )










def write_cij_energy( cij_hcp ):

    f = open('y_post_cij_energy.txt','w+')
    f.write('# Cij from energy method: \n' )

    f.write('\n%12s %12s %12s %12s %12s \n'  
        %('C11', 'C12', 'C13', 'C33', 'C44') )

    f.write('%12.4f %12.4f %12.4f %12.4f %12.4f \n' \
        %(cij_hcp[0], cij_hcp[1], cij_hcp[2], cij_hcp[3], cij_hcp[4]) )



    from myalloy import calc_elastic_constant as ce 
    E_x, E_z, nu_xy, nu_xz, mu_xz = ce.calc_transverse_isotropy( cij_hcp ) 
  

    f.write('\n%12s %12s %12s %12s %12s \n'  
        %('E_x', 'E_z', 'nu_xy', 'nu_xz', 'mu_xz') )

    f.write('%12.4f %12.4f %12.4f %12.4f %12.4f \n' \
        %( E_x,   E_z,   nu_xy,   nu_xz,   mu_xz ) )



    f.write('\n%25s %25s \n'  
        %('E_x/(1-nu_xy)', 'nu_xz/(1-nu_xy)' ) )

    f.write('%25.4f %25.4f \n' \
        %( E_x/(1-nu_xy),   nu_xz/(1-nu_xy)  ) )


    f.write('\n%25s \n'  
        %('E_x/(1-nu_xy**2)'  ) )

    f.write('%25.4f \n' \
        %( E_x/(1-nu_xy**2)   ) )




    f.close() 







main()





