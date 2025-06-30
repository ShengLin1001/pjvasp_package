import os

import shutil

import pandas as pd

from ase.io import read
from ase.io.vasp import read_vasp

import numpy as np

from scipy.interpolate import CubicHermiteSpline

from mymetal.io.general import general_read
from mymetal.io.extxyz import extxyz_to_atomlist
from mymetal.universal.plot.plot import my_plot_neb, my_plot_xy

from scipy.spatial import cKDTree
from scipy.sparse.csgraph import connected_components
from scipy.sparse import csr_matrix

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

# For detailed analysis, you can use the following function to plot the NEB results.
def analyze_neb_trajectory(file: str = None, if_save: bool = True, save_dir: str = './', cellxypoints_special_direct: list = None,
                color_list: list = [['red', 'blue', 'green', 'orange', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan',],
                        [ 'lime', 'teal', 'lavender', 'tan', 'salmon', 'gold', 'lightcoral', 'skyblue', 'darkblue', 'darkred', 'darkgreen', 'darkorange', ]],
                boundary_color: str = 'red',
                special_points_color: str = 'blue',
                fig_subp_list: list = [None, None],
                num_ax_legend_list : list = [None, None], 
                loc_list: list = ['center', 'center'],
                bbox_to_anchor_list: list = [(0.5, 0.5),(0.5, 0.5) ],
                ncol_list: list = [2, 2],
                alpha: list = [[1, 1, 0.7], [1, 1, 0.7]],
                distance_threshold: float = 0.05,
                if_plot: bool = False,):
    """Analyzes NEB (Nudged Elastic Band) trajectory data and visualizes results.

    Processes atomic trajectory data from an extxyz file, calculates displacements,
    identifies overlapping atoms, and generates visualization plots.

    Args:
        file: Path to input extxyz file.
        if_save: Whether to save results to files.
        save_dir: Directory to save output files.
        cellxypoints_special_direct: Special points in direct coordinates.
        **kwargs: Additional plotting parameters (colors, figure settings, etc.).

    Returns:
        None: Results are saved to files and optionally displayed as plots.

    Example:
        ```python
        from mymetal.post.neb import analyze_neb_trajectory
        from mymetal.io.extxyz import extxyz_to_atomlist
        from mymetal.universal.plot.general import generate_gradient_colors
        for lc in ['2.70', '2.80', '2.90']:
            workdir = f'./neb/{lc}/'
            file = workdir+"movie_CONTCAR.xyz"
            cellxypoints_special_direct = [[0, 0], [1/3, 2/3], [2/3, 1/3]]
            atomlist = extxyz_to_atomlist(file)
            num_atoms = len(atomlist[0].get_positions())
            num_frames = len(atomlist)
            color_list1 = generate_gradient_colors(if_cmap_color=True, if_reverse=True, num_colors=num_frames)
            color_list2 = generate_gradient_colors(if_cmap_color=True, cmap_color='turbo', num_colors=num_atoms)
            analyze_neb_trajectory(file, if_save=True, save_dir=workdir, cellxypoints_special_direct=cellxypoints_special_direct, color_list=[color_list1, color_list2],
                                fig_subp_list=[[5, 3], [4, 2]], num_ax_legend_list=[13, 7], ncol_list=[1,2], alpha=[[0.7, 0.7, 0.7], [0.7, 0.7, 0.7]], distance_threshold=1e-6,
                                if_plot = True)
        ```
    """

    atomlist = extxyz_to_atomlist(file)
    cell = np.array(atomlist[0].get_cell())
    cellxy = cell[:2, :2]
    cellxypoints_special = np.dot(np.array(cellxypoints_special_direct), cellxy) 
    [positions_xy_list, delta_xy_list, dist_xy_list], [positions_z_list] = get_delta_dist_list(atomlist)
    cellxypoints_xy, [figxy_width, figxy_height], axes_height_xy, figxy_lim = get_figxy_wh_lim(atomlist)
    ######################################################################################################
    # loop on every atom in every frame
    
    overlap_pairs, overlap_groups = save_overlap_pairs_groups(positions_xy_list, distance_threshold = distance_threshold, save_dir=save_dir)

    if if_plot:
        my_plot_xy(xy_list = [positions_xy_list, delta_xy_list, dist_xy_list], z_list = [positions_z_list],
                figxy_wh_lim = [cellxypoints_xy, [figxy_width, figxy_height], axes_height_xy, figxy_lim],
                fig_subp = fig_subp_list[0],
                cellxypoints_special = cellxypoints_special, if_save=if_save, save_dir=save_dir,
                color_list=color_list[0], boundary_color=boundary_color, special_points_color=special_points_color,
                num_ax_legend=num_ax_legend_list[0], loc=loc_list[0],
                bbox_to_anchor=bbox_to_anchor_list[0], ncol=ncol_list[0],
                alpha=alpha[0],)
        my_plot_xy(xy_list = [positions_xy_list, delta_xy_list, dist_xy_list], z_list = [positions_z_list],
                figxy_wh_lim = [cellxypoints_xy, [figxy_width, figxy_height], axes_height_xy, figxy_lim],
                fig_subp=fig_subp_list[1],
                cellxypoints_special = cellxypoints_special, if_save=if_save, save_dir=save_dir,
                if_every_atom=False, if_every_frame=True, if_every_frame_filename='neb_trajectory_xy_diff_frames_text.png',
                color_list=color_list[1], boundary_color=boundary_color, special_points_color=special_points_color,
                num_ax_legend=num_ax_legend_list[1], loc=loc_list[1],
                bbox_to_anchor=bbox_to_anchor_list[1], ncol=ncol_list[1],
                alpha=alpha[1],
                if_show_frame_overlap=True, overlap_groups=overlap_groups,
                )
        my_plot_xy(xy_list = [positions_xy_list, delta_xy_list, dist_xy_list], z_list = [positions_z_list],
                figxy_wh_lim = [cellxypoints_xy, [figxy_width, figxy_height], axes_height_xy, figxy_lim],
                fig_subp=fig_subp_list[1],
                cellxypoints_special = cellxypoints_special, if_save=if_save, save_dir=save_dir,
                if_every_atom=False, if_every_frame=True,
                color_list=color_list[1], boundary_color=boundary_color, special_points_color=special_points_color,
                num_ax_legend=num_ax_legend_list[1], loc=loc_list[1],
                bbox_to_anchor=bbox_to_anchor_list[1], ncol=ncol_list[1],
                alpha=alpha[1],
                )

def remove_translation(positions_list: list = None) -> list:
    """Removes global translation from a list of atomic positions.
    NOT USE!!!

    Centers each frame's positions around the center of mass of the first frame.

    Args:
        positions_list: List of atomic positions for multiple frames.

    Returns:
        list: Corrected positions with translation removed.
    """
    # center of mass
    reference_com = np.mean(positions_list[0], axis=0)  
    corrected_positions = []
    
    for positions in positions_list:
        current_com = np.mean(positions, axis=0)
        delta_com = current_com - reference_com
        corrected_positions.append(positions - delta_com)

    return corrected_positions

def get_delta_dist_list(atomlist: list = None) -> tuple:
    """Calculates displacements and distances for atoms in a trajectory.

    Args:
        atomlist: List of ASE Atoms objects representing trajectory frames.

    Returns:
        tuple: Contains (positions_xy_list, delta_xy_list, dist_xy_list) and 
               positions_z_list for all atoms across frames.
    """

    num_frames = len(atomlist)
    num_atoms = len(atomlist[0])
    cellxy = np.array(atomlist[0].get_cell())[:2,:2]
    # positions_list[j][i], jth frame, ith atom
    positions_list = []
    for j in range(num_frames):
        positions_list.append(np.array(atomlist[j].get_positions()))
    # remove the rigid translation
    # corrected_positions_list[j][i], jth frame, ith atom
    # maybe wrong
    # corrected_positions_list = remove_translation(positions_list)
    # positions_list = corrected_positions_list.copy()
    #########
    # positions_xy_list[i][j], ith atom, jth frame
    positions_xy_list = []
    positions_z_list  = []
    for i in range(num_atoms):
        positions_xy_list.append(np.array([positions_list[j][i, :2] for j in range(num_frames)]))
        positions_z_list.append(np.array([positions_list[j][i, 2] for j in range(num_frames)]))
    # 位移
    delta_xy_list = []
    # 路程
    dist_xy_list = []
    for i in range(num_atoms):
        positions_xy = np.array(positions_xy_list[i])
        delta_xy = positions_xy - positions_xy[0]
        # periodic cell 
        # A_cartesian = A_direct @ cell
        # all those martix is row-based verctor
        # A_direct in (-0.5, 0.5)
        delta_xy_direct = np.dot(delta_xy, np.linalg.inv(cellxy))
        delta_xy_direct -= np.round(delta_xy_direct)
        delta_xy = np.dot(delta_xy_direct, cellxy)
        dist_xy  = np.linalg.norm(delta_xy, axis = 1)
        delta_xy_list.append(delta_xy)
        dist_xy_list.append(dist_xy)

    return [positions_xy_list, delta_xy_list, dist_xy_list], [positions_z_list]

def get_figxy_wh_lim(atomlist: list = None) -> tuple:
    """Calculates figure dimensions and limits based on unit cell.

    Args:
        atomlist: List of ASE Atoms objects (only first frame is used).

    Returns:
        tuple: Contains (cell boundary points, figure dimensions, 
               axes height, and axis limits).
    """
    cell = np.array(atomlist[0].get_cell())
    cellxy = cell[:2, :2]
    # boundary
    cellxypoints_xy = np.array([[0, 0], cellxy[0], cellxy[0] + cellxy[1], cellxy[1], [0, 0]])
    cellxypoints_x  = cellxypoints_xy[:, 0]
    cellxypoints_y  = cellxypoints_xy[:, 1]
    figxy_width  = max(cellxypoints_x) - min(cellxypoints_x)
    figxy_height = max(cellxypoints_y) - min(cellxypoints_y)
    figxy_lim = [[min(cellxypoints_x), max(cellxypoints_x)], [min(cellxypoints_y), max(cellxypoints_y)]]
    one_fig_wh_xy = [10.72, 10.72 * figxy_height / figxy_width]
    axes_height_xy = 7.32 * figxy_height / figxy_width
    #cellz = cell[2, 2]
    return cellxypoints_xy, one_fig_wh_xy, axes_height_xy, figxy_lim

def save_overlap_pairs_groups(positions_xy_list: list = None, distance_threshold: float = 0.05, save_dir: str = './') -> tuple:
    """Identifies and saves overlapping atom pairs and groups.

    Args:
        positions_xy_list: List of xy positions for all atoms across frames.
        distance_threshold: Maximum distance to consider atoms overlapping.
        save_dir: Directory to save output files.

    Returns:
        tuple: (pairs_list, groups_list) of overlapping atoms.
    """
    pairs_list, groups_list = find_overlap_pairs_groups(positions_xy_list, distance_threshold)
    # pair, group index start from 0
    with open(os.path.join(save_dir, 'overlap_pairs.txt'), 'w') as f:
        f.write('#Pairs of overlapping atoms (index starts from 1):\n')
        for i, pairs in enumerate(pairs_list):
            f.write('\nFrame {:>4d}:\n'.format(i+1))
            for pair in pairs:
                f.write('{:>4d} {:>4d}\n'.format(pair[0]+1, pair[1]+1))

    with open(os.path.join(save_dir, 'overlap_groups.txt'), 'w') as f:
        f.write('#Groups of overlapping atoms (index starts from 1):\n')
        for i, groups in enumerate(groups_list):
            f.write('\nFrame {:>4d}:\n'.format(i+1))
            for group in groups:
                formatted_group = ' '.join('{:>4d}'.format(atom+1) for atom in group)
                f.write('  {}\n'.format(formatted_group))

    return pairs_list, groups_list
    
def find_overlap_pairs_groups(positions_xy_list: list = None, distance_threshold: float = 0.05) -> tuple:
    """Identifies overlapping atom pairs and groups using KDTree.

    Args:
        positions_xy_list: List of xy positions for all atoms across frames.
        distance_threshold: Maximum distance to consider atoms overlapping.

    Returns:
        tuple: (pairs_list, groups_list) for each frame in the trajectory.
    """
    num_frames = len(positions_xy_list[0])
    num_atoms = len(positions_xy_list)
    # pairs_list[j][i], jth frame, ith pairs
    # groups_list[j][i], jth frame, ith group of atoms
    pairs_list = []
    groups_list = []
    for i in range(num_frames):
        positions_xy = np.array([positions_xy_list[j][i] for j in range(num_atoms)])
        #distance_matrix = squareform(pdist(positions_xy))
        #print(distance_matrix)
        tree = cKDTree(positions_xy)
        pairs = tree.query_pairs(distance_threshold, output_type='ndarray')

        # 处理无重叠的情况
        if len(pairs) == 0:
            pairs_list.append( np.empty((0, 2), dtype=int))
            groups_list.append([])  # 该帧无重叠群组
            continue

        # 正常处理有重叠的情况
        pairs_list.append(pairs)
        # 找到连通区域（重叠原子群）
        # 构建邻接矩阵（稀疏矩阵）
        n_nodes = np.max(pairs) + 1  
        data = np.ones(len(pairs), dtype=bool)
        adj_matrix = csr_matrix((data, (pairs[:, 0], pairs[:, 1])), shape=(n_nodes, n_nodes))

        # 计算连通分量, labels是每个节点的连通分量标签
        # n_components是连通分量的数量
        n_components, labels = connected_components(adj_matrix, directed=False)
        # 按连通分量分组
        groups = {}
        for node, label in enumerate(labels):
            if label not in groups:
                groups[label] = []
            groups[label].append(node)
        groups_list.append([sorted(group) for group in groups.values()])
    return pairs_list, groups_list
