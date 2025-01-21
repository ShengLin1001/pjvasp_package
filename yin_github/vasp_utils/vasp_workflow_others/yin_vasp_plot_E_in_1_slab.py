#!/public1/home/a8s000114/opt/bin/python3


import yin_vasp_plot_E_in_func as ef 


ef.post_E_in(deform_type='E_in_1', struc_type='slab') 



E_in_2 = ef.get_param(filename='../y_E_in_2_slab/y_post_E_in.txt') 
E_in_1 = ef.get_param(filename='../y_E_in_1_slab/y_post_E_in.txt') *2

ef.calc_E_x_nu_xy( E_in_2, E_in_1)




