import os
from mymetal.build.bulk.planar import LIST_SSF_FAULT, get_ssf_name
from mymetal.post.planar import post_ssf

# tilted-cell 法：每胞只有 1 个层错（顶面与底面在 PBC 下是同一个面）
NFAULT = 1

nfail = 0
for lat, gsfe_type, fault, t1, t2 in LIST_SSF_FAULT:
    name = get_ssf_name(gsfe_type, fault)
    workdir = os.path.join('./y_ssf', lat, name)
    dumpdir = os.path.join(workdir, 'dump')
    print('================ 📁 %s' % workdir)
    try:
        post_ssf(post_data_path=os.path.join(dumpdir, 'y_post_data.txt'),
                 latoms_path=os.path.join(dumpdir, 'movie.lammpstrj'),
                 ideal_path=os.path.join(dumpdir, 'ideal.lammpstrj'),
                 latoms_format='lammps-dump-text',
                 nfault=NFAULT,
                 label='%s (tilted-cell, t0 = %s, %s)' % (name, t1, t2),
                 save_txt_path=os.path.join(workdir, 'y_post_ssf.txt'))
    except (SystemExit, Exception) as e:   # fail() 抛 SystemExit，缺文件抛 OSError，两类都要接住
        # 单个滑移系失败不该带垮整轮 post：其余层错照算，最后统一汇总退出码
        print('❌ %s failed: %s' % (name, e))
        nfail += 1

print('================ 📊 Summary: ssf')
print('  total %d, failed %d' % (len(LIST_SSF_FAULT), nfail))
raise SystemExit(1 if nfail else 0)
