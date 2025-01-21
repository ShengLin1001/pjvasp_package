
import os 
from myvasp import vasp_func as vf 



for file in os.listdir('./'):
    if file.startswith("bestsqs-"):
        print(file)
        vf.bestsqs_to_POSCAR(filename=file)






