

# https://en.wikipedia.org/wiki/Transverse_isotropy

import numpy as np
from myvasp import vasp_func as vf 
from myvasp import vasp_shift_to_complete_layers as vs 

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt



#==========================



def post_E_in(deform_type='E_in_2', struc_type='bulk'):

    latoms = vf.get_list_of_atoms()

    nlayers, nmiss = vs.check_layers(latoms[0])
    print( nlayers, nmiss )

    jobn, Etot, Eent, pres = vf.vasp_read_post_data()
    
    #==========================      

    natoms, a, t, tv = calc_a_t(latoms, struc_type) 

    a0, t0 = calc_a0_t0(a, Etot, t) 

    A0 = calc_A0(latoms, el=a/a0, deform_type=deform_type)
    
    if struc_type=='bulk':
        V0 = A0 * t0
    elif struc_type=='slab':
        V0 = A0 * t0 * nlayers/(nlayers-1)
    

    if deform_type == 'E_in_2':
        check_latt12_E_in_2(latoms, jobn)

    elif deform_type == 'E_in_1':
        check_latt12_E_in_1(latoms, jobn)


    
    sys_name = get_sys_name(latoms[0]) 
    
    if struc_type == 'bulk':
        check_latt3_bulk(latoms, t/(nlayers-1)*nlayers ) 
        sys_name = '%s, bulk'  %(sys_name) 

    elif struc_type == 'slab':
        print('ok1')
        try:
            check_latt3_slab(latoms)
        except Exception as e:  # Catch any exception and store the error message
            print(f"Warning: check_latt3_slab failed with error: {e}")
        sys_name = '%s, $t_v \\geqslant %.2f~\\mathrm{\AA}$' %(sys_name, tv.min() ) 
        print('ok1')
        d0 = calc_d0(latoms, elt=t/t0)
        plot_E_in_d0(sys_name, d0) 


    plot_E_in(natoms, sys_name, Etot, V0, a, a0, t, t0, deform_type)
    write_E_in(natoms, nlayers, sys_name, jobn, Etot, a, t, tv, a0, t0, A0, V0) 



#==========================





def calc_a_t(latoms, struc_type):
    natoms = latoms[0].get_positions().shape[0]

    a = np.array([])   # lattice constant 
    t = np.array([])   # slab thickness
    tv = np.array([])  # vacuum thickness 

    for i in np.arange(len(latoms)):
        latt = latoms[i].get_cell()[:] 
        a = np.append(a, latt[0,0] )


        if struc_type=='bulk':
            t  = np.append(t, latt[2,2])
            tv = np.append(tv, 0)

        elif struc_type=='slab':
            pos = latoms[i].get_positions()
            vf.confirm_0( pos.shape[0] - natoms )
            
            temp = pos[:,2].max() - pos[:,2].min()
            temp2 = latt[2,2] - temp 
            
            t = np.append(t, temp)
            tv = np.append(tv, temp2) 


    return natoms, a, t, tv     








def calc_a0_t0(a, Etot, t):
    param = np.polyfit(a, Etot, 2)
    a0 = -param[1]/(2*param[0])     # a at equilibrium 

    param2 = np.polyfit(a, t, 1)
    fun2 = np.poly1d(param2)
    t0 = fun2(a0)

    return a0, t0 







def calc_A0(latoms, el, deform_type):
    A = np.array([])
    for i in np.arange(len(latoms)):
        latt = latoms[i].get_cell()[:] 

        if deform_type == 'E_in_2':
            latt = latt / el[i]

        elif deform_type == 'E_in_1':
            latt[:,0] = latt[:,0] / el[i]

        temp = np.linalg.norm( np.cross( latt[0,:], latt[1,:] ) )
        A = np.append(A, temp) 
    
    vf.confirm_0( A.std() )
    return A.mean()





def calc_d0(latoms, elt):

    natoms = latoms[0].get_positions().shape[0]

    d = np.zeros([ len(latoms), natoms-1 ])

    for i in np.arange(len(latoms)):
        latt = latoms[i].get_positions()
        d[i,:] = np.diff(latt[:,2]) / elt[i]
               

    vf.confirm_0( np.std(d, axis=0) / np.mean(d, axis=0) *1e-8 )
    return np.mean(d, axis=0)










def get_sys_name(atoms):
    from ase.formula import Formula
    sys_name = atoms.get_chemical_formula()
    sys_name = Formula(sys_name).format('latex')
    return sys_name 






def check_latt12_E_in_2(latoms, jobn):
    latt_ref = latoms[0].get_cell()[:]    # 1st in jobn 

    for i in np.arange(len(latoms)):
        latt = latoms[i].get_cell()[:] 
        vf.confirm_0( latt[0:2,:] / latt[0,0] - latt_ref[0:2,:] / latt_ref[0,0]    )

        if jobn[i][0] == '0' or jobn[i][0] == '1' :   
            vf.confirm_0( latt[0:2,:] / (float(jobn[i])) - latt_ref[0:2,:] / (float(jobn[0]))    )






def check_latt12_E_in_1(latoms, jobn):
    latt_ref = latoms[0].get_cell()[:]    # 1st in jobn 

    for i in np.arange(len(latoms)):
        latt = latoms[i].get_cell()[:] 
        vf.confirm_0( latt[0:2, 0] / (float(jobn[i])) - latt_ref[0:2, 0] / (float(jobn[0]))    )
        vf.confirm_0( latt[0:2, 1:] - latt_ref[0:2, 1:] )






def check_latt3_bulk(latoms, t2):
 
    for i in np.arange(len(latoms)):
        latt = latoms[i].get_cell()[:] 
#        print( latt[2,:]  )
#        vf.confirm_0( latt[2,:] - np.array([0, 0, t2[i]]), 'wrong latt3 bulk' )
#        vf.confirm_0( latt[2,2] - t2[i], 'wrong latt3 bulk' )






def check_latt3_slab(latoms):
    latt_ref = latoms[0].get_cell()[:]    # 1st in jobn 

    for i in np.arange(len(latoms)):
        latt = latoms[i].get_cell()[:] 
        #vf.confirm_0( latt[2,:] - latt_ref[2,:] )




#==========================






def plot_E_in_d0(sys_name, d0):

    nd = len(d0)

    if np.mod(nd, 2) == 1:
        iref = int( (nd-1)/2 )
    else:
        iref = int( nd/2)
  

    fig_wh = [3.15, 2.5]
    fig_subp = [1, 1]
    fig1, ax1 = vf.my_plot(fig_wh, fig_subp)
     
    fig_pos  = np.array([0.26, 0.19, 0.71, 0.73])
    ax1.set_position(fig_pos)
      
    xi = np.arange(len(d0))+1

    ax1.plot(xi, d0/d0[iref]-1, '-o') 

    ax1.set_xlabel('Index $i$ of interlayer distance $d_i$')

    str1 = '$d_i / d_{ %d } -1$' %(iref+1)
    ax1.set_ylabel(str1) 

    str1 = '%s' %(sys_name)

    xlim = ax1.get_xlim()
    ylim = ax1.get_ylim()
    ax1.text( xlim[0]+(xlim[1]-xlim[0])*0.5, ylim[0]+(ylim[1]-ylim[0])*0.8, \
        str1, horizontalalignment = 'center' )


    plt.savefig('y_post_E_in.d0.pdf')
    plt.close('all')









def myfit( Etot, a, a0, V0 ):
    e = a/a0-1
     
    qe = vf.phy_const('qe')
    ed = Etot/V0 *(qe*1e21)     # energy density GPa
  
    param = np.polyfit(e, ed, 2)
    vf.confirm_0(param[1])

    fun = np.poly1d(param)
    edp = fun(e)

    R2 = calc_R2(ed, edp)

    E0 = param[-1] *V0 /(qe*1e21)     # total energy (eV)

    return e, ed, param, R2, E0  





def calc_R2(y, f):
    RES = np.sum( ( y-f        )**2 ) 
    TOT = np.sum( ( y-y.mean() )**2 ) 
    R2 = 1-RES/TOT
    return R2









def plot_E_in(natoms, sys_name, Etot, V0, a, a0, t, t0, deform_type):

    e, ed, param, R2, E0 = myfit( Etot, a, a0, V0 )

    xi = np.linspace(e[0], e[-1], 1000)
    fun = np.poly1d(param)
    yi = fun(xi)


    fig_wh = [3.15, 5]
    fig_subp = [2, 1]
    fig1, ax1 = vf.my_plot(fig_wh, fig_subp)
     
    fig_pos  = np.array([0.25, 0.55, 0.70, 0.4])
    fig_dpos = np.array([0, -0.45, 0, 0])
    for i in np.arange(2):
        ax1[i].set_position(fig_pos + i* fig_dpos)
      

    ax1[0].plot(e, ed-param[-1], 'o') 
    ax1[0].plot(xi, yi-param[-1], '-')    

    xlim = ax1[0].get_xlim()

    ylim = ax1[0].get_ylim()
    ax1[0].plot( np.array([0, 0]), ylim, '--k', alpha=0.2)
    ax1[0].set_ylim(ylim) 

    ax1[0].set_ylabel('Elastic energy density (GPa)')
    

    #==========================

    et = t/t0-1 

    param2 = np.polyfit(e, et, 1)
    vf.confirm_0(param2[-1])

    fun2 = np.poly1d(param2)
    yi2 = fun2(xi)

    ax1[1].plot(e, et, 'o') 
    ax1[1].plot(xi, yi2, '-')    


    ylim2 = ax1[1].get_ylim()
    ax1[1].plot( np.array([0, 0]), ylim2, '--k', alpha=0.2)
    ax1[1].set_ylim(ylim2) 

    ax1[1].set_ylabel('Out-of-plane strain $\\epsilon_{zz}$')
  

  

    if deform_type == 'E_in_2':

        ax1[1].set_xlabel('Biaxial in-plane strain $\\epsilon_{xx}$ and $\\epsilon_{yy}$')

        str0 = '%s\n$a_0 = %.4f~\\mathrm{\AA}$\n\n$\\frac{ E_{x} }{ 1-\\nu_{xy} } = %.2f$ GPa\n$E_0 = %.4f$ eV/atom'  \
            %(sys_name, a0, param[0], E0/natoms)

        str1 = '$t_0 = %.4f~\\mathrm{\AA}$\n\n$\\frac{ \\nu_{xz} }{ 1-\\nu_{xy} } = %.4f $' \
            %(t0, param2[0]/(-2))



    elif deform_type == 'E_in_1':
        
        ax1[1].set_xlabel('Uniaxial in-plane strain $\\epsilon_{xx}$ ($\\epsilon_{yy} = 0$)')

        str0 = '%s\n$a_0 = %.4f~\\mathrm{\AA}$\n\n$\\frac{ E_{x} }{ 1-\\nu_{xy}^2 } = %.2f$ GPa\n$E_0 = %.4f$ eV/atom'  \
            %(sys_name, a0, param[0]*2, E0/natoms)

        str1 = '$t_0 = %.4f~\\mathrm{\AA}$\n\n$\\frac{ \\nu_{xz} }{ 1-\\nu_{xy} } = %.4f $' \
            %(t0, param2[0]/(-1))


    ax1[0].text( xlim[0]+(xlim[1]-xlim[0])*0.25, ylim[0]+(ylim[1]-ylim[0])*0.55, \
        str0, horizontalalignment = 'left' )
  

    ax1[1].text( xlim[0]+(xlim[1]-xlim[0])*0.55, ylim2[0]+(ylim2[1]-ylim2[0])*0.7, \
        str1, horizontalalignment = 'left' )
  


  
    plt.savefig('y_post_E_in.pdf')
    plt.close('all')










def calc_a1(jobn, a):
    a1 = np.array([])
    for i in np.arange( len(jobn) ):
        if jobn[i][0] == '0' or jobn[i][0] == '1' :   
            a1 = np.append(a1, a[i] / (float(jobn[i])) )
    vf.confirm_0( a1.std() )
    return a1.mean()




def write_E_in(natoms, nlayers, sys_name, jobn, Etot, a, t, tv, a0, t0, A0, V0):

    e, ed, param, R2, E0 = myfit( Etot, a, a0, V0 )
    a1 = calc_a1(jobn, a)


    f = open('y_post_E_in.txt','w')
    f.write('# In-plane straining of transversely isotropic materials/slabs: \n' )

    f.write('\n%s \n'  %(sys_name) )


    f.write('\n%16s %16s %16s \n' \
        %('R2', 'a1/a0-1', 'a0 (Ang)') )

    f.write('%16.8f %16.8f %16.8f \n' \
        %( R2, a1/a0-1, a0) )



    f.write('\n%16s %16s %16s \n' \
        %('natoms', 'nlayers', 'natoms/nlayers') )

    f.write('%16.0f %16.0f %16.0f \n' \
        %( natoms, nlayers, natoms/nlayers) )

 

    f.write('\n%16s %16s %16s \n' \
        %('t0 (Ang)', 'param[0] (GPa)', 'E0 (eV/atom)') )

    f.write('%16.8f %16.8f %16.8f \n' \
        %( t0, param[0], E0/natoms) )

    f.write('\n%16s %16s \n' \
        %('A0 (Ang^2)', 'V0 (Ang^3)') )

    f.write('%16.8f %16.8f \n' \
        %( A0, V0 ) )



    f.write('\n%16s %16s %16s %16s \n' \
        %('jobn', 'Etot (eV)', 'a (Ang)', 't (Ang)') )

    for i in np.arange(len(jobn)):
        f.write('%16s %16.8f %16.8f %16.8f \n' \
            %(jobn[i], Etot[i], a[i], t[i]) )


    f.write('\n%16s %16s %16s %16s \n' \
            %('ed-ed0 (GPa)', 'a/a0-1', 't/t0-1', 'tv (Ang)' ) )

    for i in np.arange(len(jobn)):
        f.write('%16.8f %16.8f %16.8f %16.8f \n' \
            %(ed[i]-param[-1], a[i]/a0-1, t[i]/t0-1, tv[i] ) )


    f.close() 






#==========================



def get_param(filename='./y_post_E_in.txt'):
    import subprocess as sp
    mycmd = 'grep -A1 param %s | tail -1 ' %( filename )
    temp = sp.getoutput(mycmd).split()
    param = float( temp[1] ) 
    return param 
    



def calc_E_x_nu_xy(E_in_2, E_in_1):
    nu_xy = E_in_2/E_in_1-1

    E_x = E_in_2*(1-nu_xy)
    vf.confirm_0(E_in_1 - E_x/(1-nu_xy**2))

    f = open('y_post_E_in.txt','a')
    f.write('\n\n\n%16s %16s %16s %16s \n' \
        %('E_in_2', 'E_in_1', 'E_x', 'nu_xy') )
    f.write('%16.8f %16.8f %16.8f %16.8f \n' \
        %( E_in_2 ,  E_in_1 ,  E_x ,  nu_xy ) )
    f.close() 





