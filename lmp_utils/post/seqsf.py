import os
from mymetal.build.bulk.planar import LIST_SEQSF_FAULT
from mymetal.post.planar import post_ssf

# 周期堆垛序列胞每胞含 2 个等价层错（Yin 2017 §3）
NFAULT = 2

nfail = 0
for lat, fault_type, n_repeat in LIST_SEQSF_FAULT:
    workdir = os.path.join('./y_seqsf', lat, fault_type)
    dumpdir = os.path.join(workdir, 'dump')
    print('================ 📁 %s' % workdir)
    try:
        post_ssf(post_data_path=os.path.join(dumpdir, 'y_post_data.txt'),
                 latoms_path=os.path.join(dumpdir, 'movie.lammpstrj'),
                 ideal_path=os.path.join(dumpdir, 'ideal.lammpstrj'),
                 latoms_format='lammps-dump-text',
                 nfault=NFAULT,
                 label='%s (stacking sequence, n=%d)' % (fault_type, n_repeat),
                 save_txt_path=os.path.join(workdir, 'y_post_ssf.txt'))
    except (SystemExit, Exception) as e:   # fail() 抛 SystemExit，缺文件抛 OSError，两类都要接住
        print('❌ %s failed: %s' % (fault_type, e))
        nfail += 1

print('================ 📊 Summary: seqsf')
print('  total %d, failed %d' % (len(LIST_SEQSF_FAULT), nfail))
raise SystemExit(1 if nfail else 0)
