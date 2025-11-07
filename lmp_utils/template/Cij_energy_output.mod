
shell mkdir ../dump/${filename}
if "${i} == 1" then "print 'jobn   mype1(eV)   EENTRO(eV)   -stress(kB)' append ../dump/${filename}/y_post_data.txt"

# Set format is necessary to avoid scientific notation in dump file
write_dump all custom ../dump/${filename}/movie.lammpstrj id type x y z modify append yes format line "%2d %2d %16.12e %16.12e %16.12e"

variable mype1 equal c_pe1

# Obtain new stress tensor
# bar to kB
variable tmp equal pxx
variable pxx1 equal ${tmp}/1000
variable tmp equal pyy
variable pyy1 equal ${tmp}/1000
variable tmp equal pzz
variable pzz1 equal ${tmp}/1000
variable tmp equal pxy
variable pxy1 equal ${tmp}/1000
variable tmp equal pxz
variable pxz1 equal ${tmp}/1000
variable tmp equal pyz
variable pyz1 equal ${tmp}/1000

# same as VASP order
print '${jobn} $(v_mype1:%16.12e) 0.0000 in kB $(v_pxx1:%16.12e) $(v_pyy1:%16.12e) $(v_pzz1:%16.12e) $(v_pxy1:%16.12e) $(v_pyz1:%16.12e) $(v_pxz1:%16.12e)' append ../dump/${filename}/y_post_data.txt