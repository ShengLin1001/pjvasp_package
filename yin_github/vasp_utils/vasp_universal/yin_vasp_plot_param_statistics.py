
#!/home/yin/opt/bin/python3

import numpy as np 
from myvasp import vasp_func as vf
import sys



def main():
    ljobs = vasp_read_post_param()

    f = open('y_post_param_statistics.txt', 'w+')
    f.write('# VASP param statistics. OK = the same in all jobs. \n' )
    
    for key, values in  ljobs[0].items():
        
        if key == 'jobname':
            continue
        
    
        l = []
        for i in np.arange( len(ljobs) ):
            l.append( ljobs[i][key] )
            

        if key == 'EDIFF':
            temp1 = np.max(l)
        elif key == 'EDIFFG' :
            temp1 = np.min(l)
        else:
            temp1 = l[0]
    
    
        if len( np.unique(l) ) == 1:
            temp2 = 'OK'
        else:
            temp2 = 'Differ!'
    

        temp3 = ''
        if key == 'ICHARG':
            l_float = []
            for num in l:
                l_float.append(float(num))
            
            l = np.array( l_float )
            if l.max() > 10:
                temp3 = 'ICHARG > 10 ?'
    
        elif key == 'EDIFFG':
            l = np.array(l)
            mask = (l > 0)
            if  len( l[mask] ) > 0.1:
                temp3 = 'EDIFFG > 0 ?'


        f.write('\n%10s: %10s %10s    %s \n' \
            %( key, temp1, temp2, temp3) )
    
    f.close()
    





def vasp_read_post_param(filename='y_post_param.txt', \
    filename2='y_post_param_2.txt', ):
       
    ljobs = []

    f = open(filename, 'r')
    next(f)
    for line in f:
        temp = line.split()
       
        job={
            'jobname': temp[0] ,

            'ENCUT':  temp[1][6:] ,
            'ISMEAR': temp[2][7:] ,

            'EDIFF':  float( temp[3][6:] ) ,
            'ISIF':   temp[4] ,
            'EDIFFG': float( temp[5] ) ,

            'LREAL': temp[6][6:] ,

            'KP': temp[7][3:] ,

            'kmesh': temp[8],

        }
                
        ljobs.append(job)
    f.close()


    f = open(filename2, 'r')
    next(f)
    k=-1
    for line in f:
        temp = line.split()
        k += 1

        if temp[0] != ljobs[k]['jobname']:
            sys.exit('ABORT: wrong post_param_2')


        addgrid = float(temp[10]) / float(temp[9]) 
        vf.confirm_int(addgrid)
        addgrid = int(addgrid)

        ljobs[k].update({
            'ISTART':  temp[1][7:] ,
            'ICHARG':  temp[2][7:] ,
            'PREC':    temp[3][5:] ,
            'ISPIN':   temp[4][6:] ,
            'VOSKOWN': temp[5][8:] ,
            'RWIGS':   temp[6][6:] ,
            'IALGO':   temp[7][6:] ,
            'ADDGRID': addgrid ,
        })
                
    f.close()


    return ljobs

      



main()




