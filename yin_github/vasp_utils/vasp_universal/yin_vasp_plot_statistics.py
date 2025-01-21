#!/home/yin/opt/bin/python3

import numpy as np 
from myvasp import vasp_func as vf 
import sys



def main():
    jobn, Etot, Eent, pres = vf.vasp_read_post_data()
    njobs = len(jobn)
    
    if njobs < 0.5:
        sys.exit('ABORT. no jobs found. ')

    atoms = vf.get_list_of_atoms()

    E0 = np.array([])
    V = np.array([])
    V0 = np.array([])
    p = np.array([])

    for i in np.arange(njobs):
        natoms = len( atoms[i].get_positions() )

        E0 = np.append(E0, Etot[i]/natoms )

        temp = atoms[i].get_volume()
        V = np.append(V, temp )
        V0 = np.append(V0, temp/natoms  )

        p = np.append(p, np.mean( pres[i,0:3] ) )

    write_output(jobn, E0, V, V0, p)






def write_output(jobn, E0, V, V0, p):
    f = open('y_post_statistics.txt', 'w+')
    f.write('# VASP statistics of y_dir: \n' )
    
    f.write('%12s %12s %12s %12s \n' \
        %('mean', 'std', 'max-mean', 'min-mean' ) )


    a_fcc = (V0*4)**(1/3)
    a_bcc = (V0*2)**(1/3)


    data = np.vstack([E0, V0, p, a_fcc, a_bcc])
    str1 = ['E0 (eV/atom):', 'V0 (Ang^3/atom):', 'p (kBar):', 'a_fcc (Ang):', 'a_bcc (Ang):']

    for i in np.arange(data.shape[0]):
        f.write('\n %s\n' %(str1[i]) )
        f.write('%12.4f %12.4f %12.4f %12.4f \n' \
        %(data[i,:].mean(), data[i,:].std(), \
            data[i,:].max()-data[i,:].mean(), \
            data[i,:].min()-data[i,:].mean() ))




    f.write('\n\n%16s %12s %12s %12s %12s \n' \
            %('jobn', 'E0', 'V', 'V0', 'p'  ) )
    for i in np.arange( len(jobn) ):
        f.write('%16s %12.4f %12.4f %12.4f %12.4f \n' \
            %(jobn[i], E0[i], V[i], V0[i], p[i]  ) )


    f.close()
    
    


main()




