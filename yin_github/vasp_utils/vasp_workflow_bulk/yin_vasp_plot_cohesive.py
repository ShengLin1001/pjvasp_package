#!/home/yin/opt/bin/python3


import numpy as np
from ase.io.vasp import read_vasp
from myvasp import vasp_func as vf 
import math    
from scipy.interpolate import interp1d
import sys



def main():
    jobn, Etot, Eent, pres = vf.vasp_read_post_data()

    ASE_Atoms = read_vasp('../y_full_relax/CONTCAR')
    atoms_pos = ASE_Atoms.get_positions()
    natoms = atoms_pos.shape[0]
    
    V0 = np.linalg.det( ASE_Atoms.cell ) / natoms
    Etot = Etot / natoms
    

    latoms = vf.get_list_of_atoms()
    latoms2 = vf.get_list_of_outcar()

    V=np.array([]) 
    magtot = np.array([])
    for i in np.arange(len(latoms)):
        temp = latoms[i].get_volume()/natoms
        V = np.append(V, temp) 
    
        if np.abs(temp-V0) < 1E-2:#1e-10:
            print('hello')
            E0 = Etot[i]
        
        try:
            temp = latoms2[i].get_magnetic_moment()/natoms
        except:
            temp = 0
        magtot = np.append(magtot, temp)
        
    magtot = np.abs(magtot)
   
 
    # scale
    k = np.array([])
    for i in np.arange(len(jobn)):
        k = np.append(k, float(jobn[i]) )
    # if np.linalg.norm( k - (V/V0)**(1/3) ) > 1e-3:# 1e-10:
    #     sys.exit('==> ABORT. wrong scaling. ')
    

    # check
    if Etot.min() != E0:
        sys.exit('E0 is wrong. Abort!')
    
    if  Etot[-5:].std() > 5e-4 :
        print('WARNING: Eatom might be wrong!')
    
    Eatom = Etot[-1]
    Ecoh =  E0 - Eatom
   
 
    Vp, p, VB, B, p0, B0 = calc_p_B(V, Etot, V0)
    write_output(E0, Eatom, Ecoh, V0, p0, B0)
    plot_output( E0, Eatom, Ecoh, V0, k, Etot, Vp, p, VB, B, p0, B0, magtot)


    
    
def calc_p_B(V, Etot, V0):
    Vp =  V[0:-1].copy() + np.diff( V)/2
    VB = Vp[0:-1].copy() + np.diff(Vp)/2 
    
    qe = vf.phy_const('qe')
    p = -np.diff(Etot) / np.diff(V) * qe*1e21  #[GPa]
    B = -np.diff(p) / np.diff(Vp) * VB   #[GPa]
    
    fp = interp1d(Vp, p)
    fB = interp1d(VB, B)
    
    p0 = fp(V0)
    B0 = fB(V0)
    print('==> p0, B0:')
    print(p0, B0)
    
    return Vp, p, VB, B, p0, B0 
    
    
    
def write_output(E0, Eatom, Ecoh, V0, p0, B0):
    f = open("y_post_cohesive.txt", "w+")
    
    f.write("# results of cohesive energy: \n" )
    
    f.write("%16s %16s %16s \n" \
    %('E0 (eV)', 'Eatom (eV)', 'Ecoh (eV)') )
    
    f.write("%16.8f %16.8f %16.8f \n" \
    %(E0, Eatom, Ecoh) )
    
    
    f.write("\n%16s \n" \
    %('V0 (Ang^3)') )
    
    f.write("%16.6f \n" \
    %( V0 ) )
    
    
    f.write("\n%16s %16s \n" \
    %('p0 (GPa)', 'B (GPa)') )
    
    f.write("%16.6f %16.8f \n" \
    %( p0, B0 ) )
    
    f.close() 
    
    
    
def plot_output(E0, Eatom, Ecoh, V0, k, Etot, Vp, p, VB, B, p0, B0, magtot):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    
    fig_wh = [3.15, 9]
    fig_subp = [4, 1]
    fig1, ax1 = vf.my_plot(fig_wh, fig_subp)
    
    
    fig_pos  = np.array([0.20, 0.77, 0.75, 0.205])
    fig_dpos = np.array([0, -0.24, 0, 0])
    
    for i in np.arange(4):    
        ax1[i].set_position(fig_pos + fig_dpos*i)
    
    
    Elimd=math.floor( Etot.min() -0.5 ) 
    Elimu=-2*Elimd
    
    plim=-math.floor( p.min()/10) *10
    
    Blimu=math.ceil( B0*1.2 /50) *50 
    Blimd=math.floor( B.min()/50) *50 
    
    maglimu = math.ceil(magtot.max() +1e-6)
    maglimd = -0.1
    
    
    ax1[0].plot([1, 1], [Elimd, Elimu], '--k')
    ax1[0].plot(k, Etot, '-o')
    
    ax1[1].plot([0, 4], [0, 0], '--k')
    ax1[1].plot([1, 1], [-plim, plim], '--k')
    ax1[1].plot((Vp/V0)**(1/3), p, '-o')
    ax1[1].plot(1, p0, 's')
    
    ax1[2].plot([0, 4], [0, 0], '--k')
    ax1[2].plot([1, 1], [Blimd, Blimu], '--k')
    ax1[2].plot((VB/V0)**(1/3), B, '-o')
    ax1[2].plot(1, B0, 's')
    
    ax1[3].plot([1, 1], [maglimd, maglimu], '--k')
    ax1[3].plot(k, magtot, '-o')
    
    
    ax1[0].set_ylim([Elimd, Elimu])
    ax1[1].set_ylim([-plim, plim])
    ax1[2].set_ylim([Blimd, Blimu])
    ax1[3].set_ylim([maglimd, maglimu])
    
    
    plt.setp(ax1[-1], xlabel='$a/a_0$')
    plt.setp(ax1[0],  ylabel='DFT energy (eV/atom)')
    plt.setp(ax1[1],  ylabel='Pressure (GPa)')
    plt.setp(ax1[2],  ylabel='Bulk modulus (GPa)')
    plt.setp(ax1[3],  ylabel='Net magnetic moment ($\\mu_B$/atom)')
    
    
    ax1[0].text(1.5, Elimd+(Elimu-Elimd)*0.6, \
    '$E_0$ = %.4f eV \n$E_\mathrm{atom}$ = %.4f eV \n$E_\mathrm{coh}$ = %.4f eV \n\n$V_0$ = %.4f $\mathrm{\AA}^3$' \
    %(E0, Eatom, Ecoh, V0)  )
    
    
    ax1[1].text(1.5, plim*0.4, \
    'from diff: \n$p_0$ = %.1f GPa ' %(p0)  )
    
    ax1[2].text(1.5, Blimd+(Blimu-Blimd)*0.6, \
    'from diff: \n$B_0$ = %.1f GPa ' %(B0)  )
    
    plt.savefig('y_post_cohesive.pdf')

    
    
main()



