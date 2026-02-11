#!/public3/home/scg6928/mysoft/env/pyenv/dft/bin/python

from pathlib import Path
import os, shutil
from ase.io import read
import numpy as np
from mymetal.universal.plot.plot import my_plot
from mymetal.universal.plot.general import general_modify_legend
from ase.units import GPa

# 1. get y_dir directory list
# 2. use a loop to chdir every y_dir
# 3. get subdirectory list of every y_dir
# 4. post-process the subdirectory list, get force and energy, ... etc.
# 5. output the result to the root directory for every y_dir


def my_univ_post_convergence(workflow = 'y_post_convergence'):

    ly_dir = list(Path.cwd().rglob("y_dir"))
    for y_dir in ly_dir:
        os.chdir(y_dir)

        current_dir = Path.cwd()
        output_dir = os.path.join(current_dir.parent, workflow)
        try:
            shutil.rmtree(output_dir)
        except:
            print('==> no folder')

        lsub_dir = sorted(
                            [p for p in current_dir.iterdir() if p.is_dir()],
                            key=lambda p: p.name
                        )

        llenergy = []
        llforce = []
        lllstress = []

        for subdir in lsub_dir:
            # file_outcar = os.path.join(subdir, "OUTCAR")
            file_xml = os.path.join(subdir, "vasprun.xml")

            output_file = os.path.join(output_dir, f"{Path(output_dir).name}_{subdir.name}.pdf")
            os.makedirs(os.path.dirname(output_file), exist_ok=True)

            lframe, lenergy, lforce, llstress, ediffg, isif, natoms = get_infos_from_xml(file_xml)

            llenergy.append(lenergy)
            llforce.append(lforce)
            lllstress.append(llstress)

            my_plot_univ_post_convergence(lframe, lenergy, lforce, llstress, ediffg, isif, natoms, output_file)        



def get_infos_from_xml(path_xml: str = None) -> tuple:
    
    # If we need to extract more information, we can use the vasprun.xml file, which contains more detailed information about the calculation.
    # outcar = read(file_outcar, format='vasp-out', index=':')
    xml = read(path_xml, format='vasp-xml', index=':')
    nframes = len(xml)
    natoms = len(xml[0])

    # start from 1.
    lframe = list(range(1, nframes + 1))
    lenergy = [xml[iframe].get_total_energy() for iframe in range(nframes)]
    lforce = []
    llstress = []
    ediffg = xml[0].calc.parameters['ediffg']
    isif = xml[0].calc.parameters['isif']

    for iframe in range(nframes):
        force = xml[iframe].get_forces()
        # force: shape (natoms, 3)
        lforce.append(np.max(force))
        stress = xml[iframe].get_stress(voigt = True)  # eV/A^3
        # stress /= (-0.1 * GPa)                       # eV/A^3 -> kB
        lstress = stress.tolist()
        llstress.append(lstress)

    return lframe, lenergy, lforce, llstress, ediffg, isif, natoms


def my_plot_univ_post_convergence(lframe, lenergy, lforce, llstress, ediffg, isif, natoms, output_file):
    fig_subp = [3, 1]
    fig, axes = my_plot(fig_subp=fig_subp, fig_sharex = False, left = 2.6, top = 0.5)
    lax = axes.flatten()

    ax = lax[0]
    ax.set_xlabel("Frame (-)")
    ax.set_ylabel("$\Delta E$ (meV/atom)")
    lenergy_rel = [(energy - lenergy[-1]) * 1000 / natoms for energy in lenergy]
    ax.plot(lframe, lenergy_rel, '-o')

    ax = lax[1]
    ax.set_xlabel("Frame (-)")
    ax.set_ylabel("Max Force (eV/Å)")
    ax.plot(lframe, lforce, '-o')
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    ax.set_ylim(bottom=-0.01)
    ax.axhline(y=-ediffg, color='C1', linestyle='--', label=f'EDIFFG = {ediffg} eV/Å')
    general_modify_legend(ax.legend(loc='upper right',
                            bbox_to_anchor=(0.95, 0.95)),
                            alpha = 0.5) 

    ax = lax[2]
    ax.set_xlabel("Frame (-)")
    ax.set_ylabel(r"Stress (eV/Å$^3$)")
    llabel = [r'$\sigma_{xx}$', r'$\sigma_{yy}$', r'$\sigma_{zz}$', 
              r'$\sigma_{yz}$', r'$\sigma_{xz}$', r'$\sigma_{xy}$']
    lmarker = ['-o', '-s', '-^', '-d', '-v', '-<']
    for i in range(6):
        lstress_i = [lstress[i] for lstress in llstress]
        ax.plot(lframe, lstress_i, lmarker[i], label=llabel[i], alpha = 0.7)
    general_modify_legend(ax.legend(loc='upper right',
                            bbox_to_anchor=(0.95, 0.95),
                            ncol = 2),
                            alpha = 0.5)        
    ax.text(
        0.02, 0.98, f"ISIF={isif}",
        transform=ax.transAxes,
        ha="left", va="top", zorder = 10, color = 'black'
    )


    fig.savefig(output_file)

my_univ_post_convergence()