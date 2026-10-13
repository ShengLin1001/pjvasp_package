from mymetal.build.bulk.gsfe import create_gsfe_model
from mymetal.build.bulk.planar import get_ssf_size
from ase.io import write

# 与 gsfe_model.py 同一套滑移面模型；超胞尺寸按 planar.DICT_SSF_SIZE 单独取（只影响 ssf 段，
# γ-line 与表面能仍用 create_gsfe_model 的默认尺寸）。get_ssf_size 返回 None 即默认尺寸。
a = aa_template
c = cc_template

latnum = lat_template

gsfe_type = "gsfe_type_template"

atoms = create_gsfe_model(gsfe_type=gsfe_type, a=a, c=c, size=get_ssf_size(gsfe_type))

write('data.ini', atoms, format='lammps-data')
