import matplotlib.pyplot as plt

def plot_param():
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams.update({'font.size': 28})
    # cax = gamma.E_gsf_surface_plot(a1vect=[-1., -2., 0.],
    #                         a2vect=gamma.a2vect,
    #                         length_unit=length_unit,
    #                         energyperarea_unit=energyperarea_unit,
    #                         if_cbar = False)
    # fontsize = 30
    # labelsize = 28
    # # 获取当前轴对象
    # ax = plt.gca()
    # # 设置边框线宽
    # for spine in ax.spines.values():
    #     spine.set_linewidth(3)
    # plt.xlabel(r'x along <10$\overline{1}0$> ($\mathit{nm}$)', fontsize=fontsize)
    # plt.ylabel(r'y along <11$\overline{2}0$> ($\mathit{nm}$)', fontsize=fontsize)
    # # 修改刻度标签字体大小, 修改主刻度线的粗细
    # plt.tick_params(axis='both', which='major', width=2, length = 6, labelsize=labelsize, direction='in')
    # cbar = plt.colorbar(aspect=40, fraction=0.1)
    # cbar.ax.set_ylabel(f'$γ_{{gsf}}$ (${energyperarea_unit}$)', fontsize=fontsize)
    # plt.show()