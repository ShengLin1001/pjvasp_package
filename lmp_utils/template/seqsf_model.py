from mymetal.build.bulk.stacking import create_stacking_model, get_stacking_sequence
from ase.io import write

a = aa_template
c = cc_template

fault_type = "seqsf_type_template"
n_repeat = seqsf_n_template

# perfect / faulted 用同一个 (a, c) 与同一套胞约定，两者只差层序，
# 这样 E0_bulk 与面积 A 可以直接互相套用。
lat, seq_perfect, seq_fault = get_stacking_sequence(fault_type, n_repeat)
print('==> %s (%s): perfect %d layers, fault %d layers'
      % (fault_type, lat, len(seq_perfect), len(seq_fault)))

write('data.ref',   create_stacking_model(seq_perfect, a=a, c=c, lat=lat), format='lammps-data')
write('data.fault', create_stacking_model(seq_fault,   a=a, c=c, lat=lat), format='lammps-data')
