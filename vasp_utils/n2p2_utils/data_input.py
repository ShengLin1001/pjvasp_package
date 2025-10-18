import os
import StructureForge.InputTools.VASPReader as vr


# data file
data = "/store/tliu/dataset_nnp/0101/vasp/"

folders = [data + f + '/' for f in os.listdir(data)
                 if os.path.isdir(data + f)]
folders.sort()

# define output file
outf = 'input.data'

# make sure file is empty
if outf in os.listdir(os.getcwd()):
    os.remove(outf)

# collect structures
n = 0
for f in folders:
    n+=1
    fp_tmp = vr.reader(f)
    fp_tmp.write(outf, writepolicy = 'a')

print(n,'structures')


