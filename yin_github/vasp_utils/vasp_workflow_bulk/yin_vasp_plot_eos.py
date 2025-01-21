#!/home/yin/opt/bin/python3


import numpy as np
from ase.io.vasp import read_vasp
from ase.eos import EquationOfState
from myvasp import vasp_func as vf 
import sys


def main():

    latoms = vf.get_list_of_atoms()
    latoms2 = vf.get_list_of_outcar()

    natoms = latoms[0].get_positions().shape[0]

    V = np.array([]) 
    ca = np.array([])
    magtot = np.array([])

    for i in np.arange(len(latoms)):
        temp = latoms[i].get_volume()/natoms
        V = np.append(V, temp) 
        
        temp = latoms[i].get_cell()[:]
        temp = np.linalg.norm( temp[2,:] ) / np.linalg.norm( temp[0,:] ) 
        ca = np.append(ca, temp) 


        try:
            temp = latoms2[i].get_magnetic_moment()/natoms
        except:
            temp = 0
        magtot = np.append(magtot, temp)
        
    magtot = np.abs(magtot)


    jobn, Etot, Eent, pres = vf.vasp_read_post_data()

    a_pos = np.array([])  # scale in POSCAR
    for i in jobn:
        filename = './y_dir/%s/CONTCAR' %(i)
        with open(filename,'r') as f:
            a_pos = np.append( a_pos, float( f.readlines()[1] ) )

 
    Etot = Etot / natoms
    fitres = myfitting(V, Etot)


    V1 = np.array([])
    a1_pos = np.array([])
    p_dft = np.array([])

    for i in np.arange(len(jobn)):
        p_dft = np.append( p_dft, pres[i,0:2].mean()*(0.1) )  # pressure [GPa]
        
        if jobn[i][0] == '0' or jobn[i][0] == '1' :   
            V1 = np.append( V1, V[i] / float(jobn[i]))
            a1_pos = np.append( a1_pos, a_pos[i] / (float(jobn[i]))**(1/3) )
     
    
    if (V1.std() > 1e-8) or (a1_pos.std() > 1e-8) :
        print('V1, V1.std(), a1_pos.std():')
        print( V1, V1.std(), a1_pos.std()  )
        sys.exit('V1 or a1_pos is wrong. Abort!')
    else:
        V1 = V1.mean()
        a1_pos = a1_pos.mean()

    write_output(V, Etot, fitres, V1, a1_pos, p_dft)
    plot_eos(V, Etot, fitres, p_dft, V1, ca, magtot)


    ASE_eos = EquationOfState(V, Etot, eos='birchmurnaghan')
    temp = np.array([ 
    ASE_eos.fit()[1] - fitres[0], 
    ASE_eos.fit()[0] - fitres[1],
    ASE_eos.fit()[2] - fitres[2] ])
    print('==> check with ASE fitting, diff: {0} \n'.format(temp) )



#==========================

# fitting eqn
def myeqn(parameters, vol):
    E0, V0, B0, B1  = parameters

    # Birch-Murnaghan equation of state
    E= E0 + 9*V0*B0/16 *( \
    (  (V0/vol)**(2/3)-1 )**3 *B1 \
    +( (V0/vol)**(2/3)-1 )**2 *( 6- 4*(V0/vol)**(2/3) ) \
    )

    return E


def myeosp(parameters, vol):
    E0, V0, B0, B1  = parameters

    # Birch-Murnaghan equation of state
    p= 3*B0/2    \
    *( (V0/vol)**(7/3) - (V0/vol)**(5/3) )   \
    *(1 + 3/4 *(B1-4) *( (V0/vol)**(2/3) -1) )    \
    
    return p


# minimize this function
def myerrfunc(param, y, x):
    err =  y - myeqn(param, x)
    return err
   

def myfitting(x, y):
    param0 = [ y.min(), x.mean(), 1, 5]
   
    from scipy.optimize import leastsq

    plsq, cov, infodict, mesg, ier = \
    leastsq(myerrfunc, param0, args=(y, x), full_output=1 )
   
    print( '\n==> fitted parameters = {0}\n'.format(plsq) )

    ssErr = (infodict['fvec']**2).sum()
    ssTot = ((y-y.mean())**2).sum()
    R2 = 1-(ssErr/ssTot )

    fitres=np.append(plsq, R2)
    return fitres


def write_output(V, Etot, fitres, V1, a1_pos, p_dft):
    
    qe = vf.phy_const('qe')
    E0 = fitres[0]
    V0 = fitres[1]
    B0 = fitres[2] *qe*1e21
    B1 = fitres[3]
    a1_pos_new = (V0/V1)**(1/3)*a1_pos


    f = open('y_post_eos.txt','w+')
    f.write('# VASP EOS fitting: \n' )

    f.write('%16s %16s %16s %16s \n' \
    %('E0 (eV)', 'V0 (Ang^3)', 'B0 (GPa)', 'B1') )

    f.write('%16.8f %16.8f %16.8f %16.8f \n' \
    %( E0, V0, B0, B1  ) )


    f.write('\n%16s %16s %16s \n' \
    %('R2', 'a0_fcc (Ang)', 'a0_bcc (Ang)'  ) )

    f.write('%16.8f %16.8f %16.8f \n' \
    %( fitres[-1], (V0*4)**(1/3), (V0*2)**(1/3)   ) )


    f.write('\n%16s %16s %16s \n' \
    %('V1/V0-1', '(V1/V0-1)*B0', 'V1-V0' ) )

    f.write('%16.8f %16.8f %16.8f \n' \
    %( V1/V0-1, (V1/V0-1)*B0, V1-V0 ) )


    f.write('\n%16s %16s \n' \
    %('a1_pos', 'a1_pos at V0' ) )
  
    f.write("%16.8f %16.8f \n" \
    %(a1_pos, a1_pos_new ) )


    f.write('\n%16s %16s %16s %16s\n' \
    %('V (Ang^3)', 'Etot (eV)' , 'p_dft (GPa)', 'p_dft-true' ) )
  
    for i in np.arange(len(V)):
        temp = p_dft[i] - myeosp(fitres[:-1], V[i])*qe*1e21 
        f.write('%16.8f %16.8f %16.8f %16.8f\n' \
        %(V[i], Etot[i], p_dft[i], temp) )

    f.close() 



def plot_eos(V, Etot, fitres, p_dft, V1, ca, magtot):
 
    qe = vf.phy_const('qe')
    E0 = fitres[0]
    V0 = fitres[1]
    B0 = fitres[2] *qe*1e21
    B1 = fitres[3]


    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig_wh = [3.15, 9]
    fig_subp = [4, 1]
    fig1, ax1 = vf.my_plot(fig_wh, fig_subp)
     

    fig_pos  = np.array([0.23, 0.77, 0.70, 0.205])
    fig_dpos = np.array([0, -0.24, 0, 0])
    
    for i in np.arange(4):
        ax1[i].set_position(fig_pos + i* fig_dpos)
   

    xi=np.arange(V.min()/V0, V.max()/V0, 1e-4)

    ax1[0].plot(V/V0, (Etot-E0), 'o')
    yi1 = myeqn(fitres[:-1], xi*V0) -E0
    ax1[0].plot(xi, yi1, '-')

    ax1[1].plot(V/V0, p_dft, 'o')
    yi2 = myeosp(fitres[:-1], xi*V0) *qe*1e21   #[GPa]
    ax1[1].plot(xi, yi2, '-')

    ax1[2].plot(V/V0, ca, 'o')
    
    ax1[3].plot(V/V0, magtot, 'o')
    ax1[3].set_ylim([-0.1, np.ceil(magtot.max()+1e-6)])


    plt.setp(ax1[-1], xlabel='Volume $V/V_0$')
    plt.setp(ax1[0],  ylabel='Energy $E-E_0$ (eV/atom)')
    plt.setp(ax1[1],  ylabel='Pressure $p$ (GPa)')
    plt.setp(ax1[2],  ylabel='$c/a$')
    plt.setp(ax1[3],  ylabel='Net magnetic moment ($\\mu_B$/atom)')


    ax1[0].text(xi[0]+(xi[-1]-xi[0])*0.2, yi1.max()*0.7, \
    '$E_0$ = %.4f eV/atom \n$V_0$ = %.4f $\mathrm{\AA}^3$/atom \n$B_0$ = %.2f GPa \n' 
    %(E0, V0, B0)  )

    ax1[1].text(xi[0]+(xi[-1]-xi[0])*0.4, yi2.max()*0.7, \
    '$p_\mathrm{Pulay} = p_\mathrm{DFT} - p_\mathrm{true}$\n$\\approx$ %.2f GPa' 
    %( (V1/V0-1)*B0 )  )



    plt.savefig('y_post_eos.pdf')



main()



