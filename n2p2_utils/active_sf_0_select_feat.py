import numpy as np
import argparse
import os,shutil
from sklearn.preprocessing import StandardScaler
from skcosmo.feature_selection import CUR, FPS

max_r = 16
max_0 = 0.05  # max percent of 0 value SF

# select_by feat_av
def CUR_feat_av(feat_av,SFs_num,SFs_all):
    X = feat_av
    y = np.random.rand(X.shape[0])  # doesn't change result
    index = CUR(n_to_select=SFs_num).fit(X,y).selected_idx_
    SFs_set = set([SFs_all[i] for i in index])
    return SFs_set

# select_by feat_atom
#  ['Au','Ni'], 48, list of SFs_all
def CUR_feat_atom(feat_atom,ele_type,SFs_num,SFs_all):
    list_type = ele_type
    list_SFsnum = [SFs_num/len(list_type),SFs_num/len(list_type)]
    index = []
    # SF of differen ele
    for i in range(len(list_type)):
        ele = list_type[i]
        SFs_num = int(list_SFsnum[i])
        sf_index = np.where(np.array([i.split()[1] for i in SFs_all])==ele)[0]
        map_index = dict(np.c_[np.arange(len(sf_index)),sf_index])
        # selet SFs
        X = feat_atom[:,list(map_index.values())]
        # rmove other type of atom
        X = np.delete(X,np.where(X[:,0]==-1)[0],axis=0)
        y = np.random.rand(X.shape[0])  # doesn't change result
        index_ele = CUR(n_to_select=SFs_num).fit(X,y).selected_idx_
        for j in index_ele:
            index.append(map_index[j])
    SFs_set = set([SFs_all[i] for i in index])
    return SFs_set 

class SF():
    max_r = 38
#    zeta = [1,4]
#    lambd = [1]
    
    # SF('G2','Au',{'Au'}, 0.1, 3, 6)
    # SF('G3','Au',{'Au','Ni'}, 0.1, 3, 6)
    def __init__(self, G, fel, oel, eta, rs, rc, lambd=1, zeta=1):
        if(rc>self.max_r):
            raise ValueError('rc is too larger')
        if(rs>rc):
            raise ValueError('rs should small than rc')
        if(G=='G2'):
            if(len(oel)!=1):
                raise ValueError('len of oel should be 1')
        if(G=='G3'):
            if(len(oel)!=2):
                raise ValueError('len of oel should be 2')
        self.G = G
        self.fel = fel
        self.oel = oel
        self.eta = eta
        self.rs = rs
        self.rc = rc
        self.lambd = lambd
        self.zeta = zeta
        
    def _create(string):
        list_string = string.split()
        # G2
        if(list_string[2]=='2'):
            G = 'G2'
            fel = list_string[1]
            oel = [list_string[3]]
            eta = float(list_string[4])
            rs = float(list_string[5])
            rc = float(list_string[6])
            return SF( G, fel, oel, eta, rs, rc)
        
        # G3
        if(list_string[2]=='3'):
            G = 'G3'
            fel = list_string[1]
            oel = [list_string[3],list_string[4]]
            eta = float(list_string[5])
            rs = float(list_string[9])
            rc = float(list_string[8])
            zeta = float(list_string[7])
            lambd = float(list_string[6])
            return SF('G3', fel, oel, eta, rs, rc, lambd, zeta)

    def _print(self):
        if self.G=='G2':
            oel = list(self.oel)
            return 'symfunction_short %s 2 %s %.4f %.3f %.3f'%(self.fel,oel[0],self.eta,self.rs,self.rc)
        # lambda = 1, zeta = 1,4
        if self.G=='G3':
            oel = list(self.oel)
            oel.sort()
#            for l in self.lambd:
#                for z in self.zeta:
            return 'symfunction_short %s 3 %s %s %.4f %.1f %.1f %.3f %.3f'%(self.fel,oel[0],oel[1],
                                                      self.eta, self.lambd, self.zeta,self.rc, self.rs)
            
    def _explore(self):
        if self.G=='G2':
            new = []
            # eta  */ 2
            new.append(SF('G2', self.fel, self.oel, self.eta*2, self.rs, self.rc))
            new.append(SF('G2', self.fel, self.oel, self.eta/2, self.rs, self.rc))
            # rs +- 0.1*rc
            delta_rs = 0.1*self.rc
            if(self.rs+delta_rs<self.rc):
                new.append(SF('G2', self.fel, self.oel, self.eta, self.rs+delta_rs, self.rc))
            if(self.rs-delta_rs>0):
                new.append(SF('G2', self.fel, self.oel, self.eta, self.rs-delta_rs, self.rc))
            # rc */ 1.26 (2*V)
            if(self.rc*1.26<self.max_r):
                new.append(SF('G2', self.fel, self.oel, self.eta, self.rs, self.rc*1.26))
            if(self.rc/1.26>self.rs):
                new.append(SF('G2', self.fel, self.oel, self.eta, self.rs, self.rc/1.26))
                
            return new
        if self.G=='G3':
            new = []
            # eta  */ 2
            new.append(SF('G3', self.fel, self.oel, self.eta*2, self.rs, self.rc))
            new.append(SF('G3', self.fel, self.oel, self.eta/2, self.rs, self.rc))
            # rs +- 0.1*rc
            delta_rs = 0.1*self.rc
            if(self.rs+delta_rs<self.rc):
                new.append(SF('G3', self.fel, self.oel, self.eta, self.rs+delta_rs, self.rc))
            if(self.rs-delta_rs>0):
                new.append(SF('G3', self.fel, self.oel, self.eta, self.rs-delta_rs, self.rc))
            # rc */ 1.26 (2*V)
            if(self.rc*1.26<self.max_r):
                new.append(SF('G3', self.fel, self.oel, self.eta, self.rs, self.rc*1.26))
            if(self.rc/1.25>self.rs):
                new.append(SF('G3', self.fel, self.oel, self.eta, self.rs, self.rc/1.26))
            
            return new







# feature for frm
feat_atom = np.load('feat_atom.npy')
feat_av = np.load('feat_av.npy')

# remove SFs which is always zero
del_c = []
for i in range(feat_atom.shape[1]):
#    if np.array([i==0 or i==-1 for i in feat_atom[:,i]]).prod():  # every atom is 0
#    for j in range(len(feat_atom[:,i])):
#        if(feat_atom[:,i][j] != 0 and feat_atom[:,i][j] != -1):  
#            break;  
#        if j == (len(feat_atom[:,i])-1):
    if(len(np.where(feat_atom[:,i]==0)[0])/len(np.where(feat_atom[:,i]!=-1)[0])>max_0):
            del_c.append(i)
# no 0 columns
feat_atom = np.delete(feat_atom, del_c, axis=1)                    # del this column from feat_atom, feat_av
feat_av = np.delete(feat_av, del_c, axis=1)                    # del this column from feat_atom, feat_av
np.save('feat_atom.npy',feat_atom)
np.save('feat_av.npy',feat_av)

# write SFs and droped sf
if len(del_c) != 0 :
    # del from SFs_all.dat
    # every SFs
    c = 0
    # drop out SFs
    f0 = open('SFs_drop.dat','w')
    with open('SFs_all.dat','r') as f:
        lines = f.readlines() 
    with open('SFs_all.dat','w') as f:
        for line in lines:
            if(c not in del_c):
                f.write(line)
            else:
                f0.write(line)
            c += 1
    f0.close()
print('delete: ',len(del_c),' features')

# read SFs
SFs_all = []
with open('SFs_all.dat','r') as f:
    for line in f:
        SFs_all.append(line)


# num of SFs  to select
ap = argparse.ArgumentParser()
ap.add_argument("number", nargs='?', default="36")
a = ap.parse_args()
SFs_num = int(a.number)


SFs_select = list(CUR_feat_av(feat_av,SFs_num,SFs_all))
#SFs_select = list(CUR_feat_atom(feat_atom,['Au','Ni'], SFs_num, SFs_all))
SFs_select.sort()

# input.nn without SFs
with open('input.nn.all','r') as r:
    lines = r.readlines()
with open('input.nn','w') as w:
    for l in lines:
        if l[0] == '#':
            w.write(l)
        elif 'symfunction_short' not in l:
            w.write(l)
        else:
            pass
#append symmetry function to input.nn
with open('input.nn','a+') as w:
    for l in SFs_select:
        w.write(l)


# active SF
#SFs_all_set = set(SFs_all)
#SFs_new = set()
#for sf in SFs_select:
    # create a SF from string
#    sf_0 = SF._create(sf)
#    sf_0.max_r = max_r 
#    new_list = sf_0._explore()
#    [SFs_new.add(i._print()) for i in new_list]

#print('explore %d'%len(SFs_new))
# delete repeated(both itself and SFs_all) SFs
#SFs_new = SFs_new - SFs_all_set
#print('finally explored %d'%len(SFs_new))

# cp active_sf_0 to active_sf_00
#current_dir = os.getcwd().split('/')[-1]
#new_dir = current_dir + '0'
#os.chdir('../')
#os.system('cp -r %s %s'%(current_dir,new_dir))
#os.chdir('%s'%new_dir)

# earse SFs in input.nn.all
#with open('input.nn.all','r') as r:
#    lines = r.readlines()
#with open('input.nn.all','w') as w:
#    for l in lines:
#        if l[0] == '#':
#            w.write(l)                                                                                                
#        elif 'symfunction_short' not in l:
#            w.write(l)
#        else:
#            pass
# write SFs_new to input.nn.all
#with open('input.nn.all','a') as w:
#    for l in list(SFs_new):
#        w.write(l)
#        w.write('\n')


