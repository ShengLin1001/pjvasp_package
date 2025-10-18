import numpy as np
from yaya_quant.re import io
import os


# every SF
file_list = ['work_dir/'+f for f in os.listdir('work_dir/')]
file_list.sort()

# feat_atom  feat_frm
feat_atom = []
feat_av = []
# read SFs   elements, value
for f in file_list:
    os.chdir(f)
    print(f)
    fname = 'function.data'
    with open(fname,'r') as f:
        counter = 0
        i_frm = 0
        sp = []
        feat = []
        sp_frm = []
        feat_frm = []
        prop_frm = []
        nat_frm = []
        for line in f:                                                   
            split = line.split()
		# number of atoms in this frame
            if len(split) == 1:
                nat = int(split[0])
                nat_frm.append(nat)
		# counting
            elif counter < nat:
                sp.append(int(split[0]))  					# kind of atom
                feat.append(np.array([float(s) for s in split[1:]]))		# collect features
                counter += 1
		# last line of every frame
            else:
                i_frm += 1                                                     # number of structures
                sp = np.array(sp)
                feat = np.array(feat,dtype=object)
                sp_frm.append(sp)
                feat_frm.append(feat)
                sp = []
                feat = []
                counter = 0
        
        n_frm = i_frm
        sp_frm = np.array(sp_frm,dtype=object)
        feat_frm = np.array(feat_frm,dtype=object)
        nat_frm = np.array(nat_frm,dtype=object)

        for i_frm in range(n_frm):
            for i_atom in range(len(sp_frm[i_frm])):
                # 为了在多元素体系中提取出施加了单个对称函数和一个g3的原子
                if len(feat_frm[i_frm][i_atom])==2:        # 为什么要寻找只有两个对称函数的原子？对称函数的个数不是固定的吗，对每个原子都一样
                    find_sp = sp_frm[i_frm][i_atom]        # 因为他将每个对称函数分别测试，并添加了一个g3，防止n2p2报错 （in sub_cal.py）
                    break 
# column of this two matrix
        feat_atom_c = []
        feat_av_c = [] 
        for i_frm in range(n_frm):
# number of sp in this frm
            SF_sum = 0
            for i_atom in range(len(sp_frm[i_frm])):
                if sp_frm[i_frm][i_atom] == find_sp:
                    # SF: single structure, single atom, single desired test SF (max to exclude uselesss added g3)
                    SF = max(feat_frm[i_frm][i_atom])  
                    # 0 0 is 0   0 0.1 is 0.1
                    feat_atom_c.append(SF)                  # SF value of this ith frame, this jth atom
                    SF_sum += SF
                else:
                    feat_atom_c.append(-1)
            feat_av_c.append(SF_sum/len(sp_frm[i_frm]))     # SF average value of this ith frame
    feat_atom.append(feat_atom_c)   # [] every frame, every find_sp atom, single SF 
    feat_av.append(feat_av_c)       # [] every frame, average value of single SF
    os.chdir('../../') 
           
#


#io.save_pkl(feat_atom,'feat_atom.pkl')
#io.save_pkl(feat_av,'feat_av.pkl')
feat_atom = np.array(feat_atom).T
feat_av = np.array(feat_av).T
if('feat_atom.npy' not in os.listdir()):
    np.save('feat_atom.npy', feat_atom)
    np.save('feat_av.npy', feat_av) 
else:
    feat_atom_pre = np.load('feat_atom.npy')
    feat_av_pre = np.load('feat_av.npy') 
    feat_atom = np.c_[feat_atom_pre,feat_atom] 
    feat_av = np.c_[feat_av_pre,feat_av]
    np.save('feat_atom.npy', feat_atom)
    np.save('feat_av.npy', feat_av) 

# clear workdir
#os.system('rm -r work_dir/*')
   
