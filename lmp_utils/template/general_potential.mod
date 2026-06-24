# NOTE: This script can be modified for different pair styles 
# See in.elastic for more info.

# Set per-type masses (由 model.py 的 self.lele 生成；pei_lmp_run_properties 在各 workdir 写入 general_mass.mod)
# 放在 pair 之前：general_potential.mod 总是在建/读结构之后被 include，此时原子类型已定义，mass 合法。
include general_mass.mod

# Choose potential
pair_style pair_style_template
pair_coeff * * pair_coeff_template

# Setup neighbor style
neighbor        2.0 bin
neigh_modify  every  1  delay  0  check yes

# potential energy
compute         peatom1 all pe/atom
compute         pe1     all reduce   sum c_peatom1

# entropy need install EXTRA-COMPUTE package
#compute         enatom1 all entropy/atom
#compute         en1     all reduce   sum c_enatom1

# Setup minimization style
#  min_style	     cg
#  min_modify	     dmax ${dmax} line quadratic

# Setup output
thermo         100
thermo_style    custom step pe etotal press pxx pyy pzz pxy pxz pyz vol c_pe1
thermo_modify   format line "%5d  %15.10f %15.10f %15.10f %13.4e  %13.4e %13.4e %13.4e  %13.4e %13.4e %15.10f %15.10f"


