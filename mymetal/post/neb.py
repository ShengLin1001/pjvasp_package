import os
import shutil
import pandas as pd
from ase.io.vasp import read_vasp
from mymetal.universal.plot.plot import my_plot_neb

import numpy as np
from scipy.interpolate import CubicHermiteSpline
from mymetal.io.general import general_read

def post_neb(dirsurf: str = 'y_neb', refcontcar: str='./y_full_relax_neb/ini/CONTCAR', repost: bool = True):
    """Processes NEB calculation results and generates analysis files.

    Args:
        dirsurf: Directory containing NEB results. Defaults to 'y_neb'.
        refcontcar: Path to reference CONTCAR file. Defaults to ini/CONTCAR.
        repost: If True, reprocesses all results. Defaults to True.

    Returns:
        None: Generates output files and plots in the working directory.
    """
    # root directory
    myroot = os.getcwd()
    refcontcar = os.path.join(myroot, refcontcar)
    dirsurf    = os.path.join(myroot, dirsurf)
    workdir    = os.path.join(dirsurf, 'y_dir')

    # check if the directory exists
    if os.path.isdir(dirsurf):
        print(f"✅ Directory {dirsurf} exists.")
    else:
        raise FileNotFoundError(f"❌ Directory {dirsurf} does not exist. Please run the convergence calculations first.")

    # ref atoms, natoms could be read in OUTCAR (NIONS)
    # import subprocess
    # output = subprocess.check_output("grep NIONS OUTCAR | tail -1 | awk '{print $NF}'", shell=True, text=True).strip()
    # nions = int(output)
    # print(nions)
    atoms_ref = read_vasp(refcontcar)
    natoms    = atoms_ref.get_positions().shape[0]

    # general post and move / change the file
    os.chdir(dirsurf)
    if repost:
        os.system("(pei_vasp_univ_clean_all_files && rm -rf vaspgr) > /dev/null 2>&1")
        # general post
        # you need to module load gunplot
        os.chdir(workdir)
        os.system("(module load gnuplot/5.2.8 && pei_vasp_univ_neb_nebresults.pl && pei_vasp_univ_neb_nebmovie 0 && pei_vasp_univ_neb_nebmovie 1)")
        os.system("mv ./*.dat ./*.eps ./vaspgr ./movie* ../")
        os.chdir(dirsurf)

    my_copy_neb_files()

    my_add_head()

    nebdf, nebefdf, spline_df, extsdf = my_read_neb()

    mysplinedf = my_spline_neb(nebdf, nebefdf)

    my_plot_neb(nebdf, nebefdf, spline_df, extsdf, mysplinedf, natoms, if_save=True)

    my_write_neb(mysplinedf)

def my_copy_neb_files() -> None:
    """Copies and converts NEB result files to standardized formats.
    
    Converts EPS to PDF and adds prefix to DAT files.
    """
    for fname in os.listdir('vaspgr'):
        if fname.endswith(".eps"):
            eps_path = os.path.join('vaspgr', fname)
            pdf_path = os.path.splitext(eps_path)[0] + ".pdf"
            os.system(f"ps2pdf {eps_path} {pdf_path}")

    for fname in os.listdir('./'):
        if fname.endswith(".eps"):
            eps_path = os.path.join('.', fname)
            pdf_path = os.path.splitext(eps_path)[0] + ".pdf"
            os.system(f"ps2pdf {eps_path} {pdf_path}")

    for filename in os.listdir('./'):
        if filename.endswith(".dat") and not filename.startswith("p_"):
            new_name = "p_" + filename
            src_path = os.path.join('./', filename)
            dst_path = os.path.join('./', new_name)
            shutil.copyfile(src_path, dst_path)

def my_write_neb(mysplinedf: pd.DataFrame = None, outname: str = 'p_myspline.txt') -> None:
    """Writes spline interpolation results to formatted text file.

    Args:
        mysplinedf: DataFrame containing spline data.
        outname: Output filename. Defaults to 'p_myspline.txt'.
    """
    f = open(outname,'w')
    coord = mysplinedf['dist_cum(Å)']
    energy = mysplinedf['energy(eV)']
    force = mysplinedf['force_b(eV/Å)']
    f.write('# my spline data: \n' )

    f.write('\n%16s %16s %16s\n' \
        %('dist_cum(Å)', 'energy(eV)', 'force_b(eV/Å)') )

    for i in np.arange(len(coord)):
        f.write('%16.8f %16.8f %16.8f\n' \
            %(coord[i], energy[i], force[i]) )

    f.close()  

def my_add_head(files: list = ['p_neb.dat', 'p_nebef.dat', 'p_spline.dat', 'p_exts.dat'],
                heads = {
                'p_neb.dat': "{:>3} {:>16} {:>16} {:>16} {:>6}".format("id", "dist_cum(Å)", "energy_rel(eV)", "force_b(eV/Å)", "image"),
                'p_nebef.dat': "{:>4} {:>16} {:>16} {:>16}".format("id", "force_max(eV/Å)", "energy(eV)", "energy_rel(eV)"),
                'p_spline.dat': "{:>4}{:>15}{:>20}{:>20}".format("id", "dist_cum(Å)", "energy_rel(eV)", "force_b(eV/Å)"),
                'p_exts.dat': "{:8}{:2}{:>6}{:>3}{:>6}{:>10}{:>5}{:>8}{:>16}".format("t1", "id", "t2", "t3", "t4", "image", "t5", "t6", "energy_rel(eV)")
                }
                ) -> None:
    """Adds formatted headers to NEB data files.

    Args:
        files: List of files to process.
        heads: Dictionary mapping filenames to header strings.
    """
    for file in files:
        if os.path.isfile(file):
            header = heads[file]
            os.system(f"echo '{header}' | cat - {file} > tmp && mv tmp {file}")
        else:
            print(f"File {file} does not exist. Skipping header addition.")

def my_read_neb(nebfile: str = 'p_neb.dat', nebeffile: str = 'p_nebef.dat', 
                splinefile: str = 'p_spline.dat', extsfile: str = 'p_exts.dat') -> tuple:
    """Reads NEB result files into DataFrames.

    Args:
        nebfile: Path to neb.dat file.
        nebeffile: Path to nebef.dat file.
        splinefile: Path to spline.dat file.
        extsfile: Path to exts.dat file.

    Returns:
        tuple: (nebdf, nebefdf, spline_df, extsdf) DataFrames.
    """
    nebdf = general_read(nebfile, has_header=True)
    nebefdf = general_read(nebeffile, has_header=True)
    spline_df = general_read(splinefile, has_header=True)
    extsdf = general_read(extsfile, has_header=True)
    return nebdf, nebefdf, spline_df, extsdf

def my_spline_neb(nebdf: pd.DataFrame = None, nebefdf: pd.DataFrame =  None,) -> pd.DataFrame:
    """Generates cubic spline interpolation of NEB path.

    Args:
        nebdf: DataFrame containing path coordinates.
        nebefdf: DataFrame containing energy values.

    Returns:
        pd.DataFrame: Interpolated path with 1000 points.
    """
    # read data
    coord = nebdf['dist_cum(Å)']
    energy = nebefdf['energy(eV)'] 
    force = nebdf['force_b(eV/Å)']

    # my spline
    my_spline = CubicHermiteSpline(coord, energy, -force)
    my_spline_coord = np.linspace(coord.iloc[0], coord.iloc[-1], 1000)
    my_spline_energy = my_spline(my_spline_coord)
    my_spline_force = -my_spline.derivative()(my_spline_coord)

    my_spline_df = pd.DataFrame({
        'dist_cum(Å)': my_spline_coord,
        'energy(eV)': my_spline_energy,  
        'force_b(eV/Å)': my_spline_force
    })
    return my_spline_df
