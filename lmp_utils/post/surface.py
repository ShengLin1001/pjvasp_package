import os
from mymetal.build.bulk.planar import LIST_SURF_TYPE
from mymetal.post.planar import post_surface

nfail = 0
for lat, gsfe_type in LIST_SURF_TYPE:
    workdir = os.path.join('./y_surface', lat, gsfe_type)
    dumpdir = os.path.join(workdir, 'dump')
    print('================ 📁 %s' % workdir)
    try:
        post_surface(post_data_path=os.path.join(dumpdir, 'y_post_data.txt'),
                     latoms_path=os.path.join(dumpdir, 'movie.lammpstrj'),
                     latoms_format='lammps-dump-text',
                     label=gsfe_type,
                     save_txt_path=os.path.join(workdir, 'y_post_surface.txt'))
    except (SystemExit, Exception) as e:   # fail() 抛 SystemExit，缺文件抛 OSError，两类都要接住
        print('❌ %s failed: %s' % (gsfe_type, e))
        nfail += 1

print('================ 📊 Summary: surface')
print('  total %d, failed %d' % (len(LIST_SURF_TYPE), nfail))
raise SystemExit(1 if nfail else 0)
