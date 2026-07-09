"""
dataset module

This module provides functions for creating, reading, and writing datasets in NNP format.
It includes functions to extract atomic configurations, energies, and forces from VASP OUTCAR files,
as well as to write these configurations into a dataset file compatible with NNP tools.

Classes:
    - nnpdata: A class for handling NNP input datasets.

Functions:
    - load_from_datafile(file_path): Load dataset from an existing NNP input data file
    - load_from_outcar(outcarfile, index, tag, comment_file): Load dataset from a VASP OUTCAR file
    - write(outfile_name, append): Write the dataset to an NNP input data file
    - get_dict(): Get the internal dataset as a dictionary
    - read_dft_reference(dir_dft_root, ...): Batch-read DFT (VASP) reference
      stretch/cij/gsfe properties into one dict for comparison with NNP epoch scans.

Change Log:
    - 2025.10.19 Integrated all functions into the nnpdata class for better encapsulation.
    - 2026.06.23 Added read_dft_reference() to collect DFT property baselines.
"""


from ase.calculators.singlepoint import SinglePointCalculator
from ase.io.vasp import read_vasp_out
from ase import Atoms
import os
import numpy as np
import pandas as pd
import re
from collections import defaultdict
from pathlib import Path


def _choose_inplane_rep(n0: int, size_top: int = 72, size_bottom: int = 36,
                        size_close: int = 48) -> tuple:
    """Pick an in-plane (a, b) supercell factor to normalize the atom count.

    Returns ``(nx, ny)`` so that ``n0 * nx * ny`` lies in
    ``[size_bottom, size_top]`` and is as close as possible to ``size_close``.
    Replication can only grow the atom count, so structures already larger than
    ``size_top`` are kept untouched (returns ``(1, 1)``).

    Args:
        n0 (int): Atom count of the original structure.
        size_top (int): Upper bound of the allowed atom count.
        size_bottom (int): Lower bound of the allowed atom count.
        size_close (int): Target atom count to approach.

    Returns:
        tuple: ``(nx, ny)`` in-plane replication factors (``nz`` stays 1).

    Notes:
        For any ``n0 <= size_top`` a valid factor always exists: when
        ``n0`` already falls in ``[size_bottom, size_top]`` then ``k=1`` is a
        candidate, and when ``n0 < size_bottom`` the window length
        ``size_top - size_bottom`` is wide enough (>= n0) to contain a multiple.
        The chosen total multiplier ``k`` is then split into the most balanced
        ``nx * ny`` (a prime ``k`` degenerates to a ``1 x k`` in-plane strip).
    """
    if n0 > size_top:
        return (1, 1)

    best_k = None
    best_diff = None
    for k in range(1, size_top // n0 + 1):
        n = n0 * k
        if n < size_bottom or n > size_top:
            continue
        diff = abs(n - size_close)
        if best_diff is None or diff < best_diff:
            best_diff = diff
            best_k = k

    if best_k is None:  # 理论上 n0 <= size_top 时不会发生，留作兜底
        return (1, 1)

    # 把 k 拆成尽量均衡的 nx * ny（nx <= ny）
    nx = 1
    for d in range(int(best_k ** 0.5), 0, -1):
        if best_k % d == 0:
            nx = d
            break
    ny = best_k // nx
    return (nx, ny)


class nnpdata:
    """A class for handling NNP input datasets.

    This class supports three initialization methods:
      1. (Full) Provide a file path to load an existing NNP dataset.
      2. (Part) Read from an OUTCAR file and generate the dataset.

    Typically, multiple partial datasets can be written incrementally using
    the second method (`load_from_outcar`), and later loaded all at once using
    the first method (`load_from_datafile`).

    Methods:
        - load_from_datafile(file_path): Load dataset from an existing NNP input data file
        - load_from_outcar(outcarfile, index, tag, comment_file): Load dataset from a VASP OUTCAR file
        - write(outfile_name, append): Write the dataset to an NNP input data file
        - get_dict(): Get the internal dataset as a dictionary

    Example:
        ```python
        inputdata1 = nnpdata()
        inputdata1.load_from_outcar(outcarfile='./workdir2/OUTCAR', index=':', tag='test', comment_file=None)
        inputdata1.write(outfile_name='./input-test3.data')
        mydict = inputdata1.get_dict()
        ```

    Attributes:
        - latoms: List of ASE Atoms objects representing atomic configurations.
        - ltags: List of tags corresponding to each configuration.
        - lcomment_files: List of comment files for each configuration.
        - lfiles: List of source files for each configuration.
        - lstruct_numbers: List of structure numbers for each configuration.
        - lfull_struct_numbers: List of total structure counts for each configuration.
        - lforces: List of forces arrays for each configuration.
        - lenergies: List of energies for each configuration.
        
    Todo:
        - Add summary/statistics utilities for dataset inspection.
    """
    # 重复 load 会覆盖之前的数据
    def __init__(self):
        self.latoms = []
        self.ltags = []
        self.lcomment_files = []
        self.lfiles = [] # source files
        self.lstruct_numbers = []
        self.lfull_struct_numbers = []
        self.lforces = []
        self.lenergies = []

        # Full content dict
        self.dict = {}


    # n2p2 to dict (include the atoms object)
    # Full dataset read from input.data
    def load_from_datafile(self, file_path: str= None):
        """Reads a full NNP dataset from an existing input.data file.

        Args:
            file_path (str): Path to the dataset file (e.g., 'input.data').

        Raises:
            ValueError: If `file_path` is None.
            ValueError: If atomic positions, forces, and symbols have mismatched lengths.

        Example:
            >>> data = nnpdata()
            >>> data.load_from_datafile("./input.data")
        """
        if file_path is None:
            raise ValueError("The input file path is None.")
        
        with open(file_path, 'r') as f:
            lines = f.readlines()
        i = 0
        n = len(lines)

        latoms = []
        ltags  = []
        lcomment_files = []
        lfiles = [] # source files
        lstruct_numbers = []
        lfull_struct_numbers = []
        lforces = []
        lenergies = []

        while i < n:
            line = lines[i].strip()
            if line.startswith("begin"):
                lattice = []
                positions = []
                forces = []
                symbols = []
                energy = 0.0
                tag = 'all'
                source_file = 'unknown'
                struct_number = 1
                full_struct_number = 1
                begin_i = i
                i = i + 1
                while i < n and lines[i].strip() != "end":
                    line = lines[i].strip()
                    if line.startswith('comment'):
                        # Example format: comment | tag=xx | frame=1/12 | file=xxx
                        comment = line
                        tag_match = re.search(r'tag=([^\s]+)', comment)
                        struct_match = re.search(r'frame=(\d+)(?:/(\d+))?', comment)
                        file_match = re.search(r'file=([^\s]+)', comment)
                        tag = tag_match.group(1) if tag_match else 'all'
                        struct_number = int(struct_match.group(1)) if struct_match else 1
                        full_struct_number = int(struct_match.group(2)) if struct_match and struct_match.group(2) else 1
                        source_file = file_match.group(1) if file_match else 'unknown'
                        # print(tag, struct_number, full_struct_number, source_file)
                    elif line.startswith("lattice"):
                        # 1,2,3
                        lattice.append([float(x) for x in line.split()[1:4]])
                    elif line.startswith("atom"):
                        positions.append([float(x) for x in line.split()[1:4]])
                        symbols.append(line.split()[4])
                        forces.append([float(x) for x in line.split()[7:10]])
                    elif line.startswith("energy"):
                        energy = float(line.split()[1])
                    elif line.startswith("charge"):
                        charge = float(line.split()[1])
                    i = i + 1
                if not (len(positions) == len(forces) == len(symbols)):
                    raise ValueError(f"Error: positions, forces, and symbols lengths do not match in index {begin_i}.")
                atoms = Atoms(
                    symbols=symbols,
                    positions=np.array(positions),
                    cell=np.array(lattice),
                    pbc=True,
                )
                # Used to remember the energy, force and stress for a given configuration.
                atoms.calc = SinglePointCalculator(atoms, energy=energy, forces=np.array(forces), charges=charge) 
                latoms.append(atoms)
                ltags.append(tag)
                lcomment_files.append(None)
                lfiles.append(source_file)
                lstruct_numbers.append(struct_number)
                lfull_struct_numbers.append(full_struct_number)
                lforces.append(np.array(atoms.get_forces(apply_constraint=False)))
                lenergies.append(atoms.get_potential_energy())
            i = i+1
        
        self.latoms = latoms
        self.ltags = ltags
        self.lcomment_files = lcomment_files
        self.lfiles = lfiles
        self.lstruct_numbers = lstruct_numbers
        self.lfull_struct_numbers = lfull_struct_numbers
        self.lforces = lforces
        self.lenergies = lenergies

        self.dict = {
            'latoms': latoms,
            'ltags': ltags,
            'lcomment_files': lcomment_files,
            'lfiles': lfiles,
            'lstruct_numbers': lstruct_numbers,
            'lfull_struct_numbers': lfull_struct_numbers,
            'lforces': lforces,
            'lenergies': lenergies
        }

    # For one OUTCAR file, load dataset
    def load_from_outcar(self, outcarfile: str=None, index: str = ':', tag: str = 'all', comment_file: str = None):
        """Reads structures, energies, and forces from a VASP OUTCAR file.

        Args:
            outcarfile (str): Path to the OUTCAR file.
            index (str, optional): Frames to read. Default ':' means all frames.
            tag (str, optional): Tag label for the dataset, e.g., 'train' or 'test'.
            comment_file (str, optional): Optional comment label or external note file.

        Notes:
            The read data are stored internally in lists (`latoms`, `lforces`, etc.)
            and can later be written out using `write()`.

        Example:
            >>> data = nnpdata()
            >>> data.load_from_outcar("OUTCAR", index=":", tag="Mg_bulk")
        """
        # Read whole atoms configurations from OUTCAR
        # index=':' means read all configurations
        outcar = read_vasp_out(outcarfile, index=index)
        num = len(outcar)

        latoms = []
        ltags  = []
        lcomment_files = []
        lfiles = [] # source files
        lstruct_numbers = []
        lfull_struct_numbers = []
        lforces = []
        lenergies = []

        for i, atoms in enumerate(outcar):
            latoms.append(atoms)
            ltags.append(tag)
            lcomment_files.append(comment_file)
            lfiles.append(outcarfile)
            lstruct_numbers.append(i+1)
            lfull_struct_numbers.append(num)
            lenergies.append(atoms.get_potential_energy())
            lforces.append(atoms.get_forces(apply_constraint=False))

        self.latoms = latoms
        self.ltags = ltags
        self.lcomment_files = lcomment_files
        self.lfiles = lfiles
        self.lstruct_numbers = lstruct_numbers
        self.lfull_struct_numbers = lfull_struct_numbers
        self.lforces = lforces
        self.lenergies = lenergies

        self.dict = {
            'latoms': latoms,
            'ltags': ltags,
            'lcomment_files': lcomment_files,
            'lfiles': lfiles,
            'lstruct_numbers': lstruct_numbers,
            'lfull_struct_numbers': lfull_struct_numbers,
            'lforces': lforces,
            'lenergies': lenergies
        }

    def adjust_size(self, size_top: int = 72, size_bottom: int = 36, size_close: int = 48):
        """In-place: replicate each loaded structure in-plane to normalize atom count.

        For every Atoms in ``self.latoms`` an in-plane ``(nx, ny)`` factor is
        chosen via :func:`_choose_inplane_rep` so that the atom count lands in
        ``[size_bottom, size_top]`` and approaches ``size_close``; ``nz`` stays 1
        to protect the vacuum direction of slabs. Structures already larger than
        ``size_top`` are left untouched.

        On replication (``M = nx * ny``) the total energy scales as ``E * M`` and
        the forces are tiled ``np.tile(F, (M, 1))`` (matching ASE ``repeat``'s
        block-wise atom order); a fresh ``SinglePointCalculator`` is attached
        because ``repeat`` drops the calculator. ``lforces`` / ``lenergies`` and
        ``self.dict`` are updated in sync so ``write()`` and the npy dumps agree.

        Args:
            size_top (int): Upper bound of the allowed atom count.
            size_bottom (int): Lower bound of the allowed atom count.
            size_close (int): Target atom count to approach.
        """
        for i, atoms in enumerate(self.latoms):
            n0 = len(atoms)
            nx, ny = _choose_inplane_rep(n0, size_top=size_top,
                                         size_bottom=size_bottom, size_close=size_close)
            M = nx * ny
            if M == 1:
                continue

            energy = atoms.get_potential_energy()
            forces = atoms.get_forces(apply_constraint=False)

            sc = atoms.repeat((nx, ny, 1))
            new_forces = np.tile(forces, (M, 1))
            new_energy = energy * M
            sc.calc = SinglePointCalculator(sc, energy=new_energy, forces=new_forces)

            self.latoms[i] = sc
            self.lforces[i] = new_forces
            self.lenergies[i] = new_energy
            tag = self.ltags[i] if i < len(self.ltags) else 'all'
            print(f"  ↗ resize {n0}→{n0 * M} (×{nx}×{ny}) [{tag}]")

        self.dict = {
            'latoms': self.latoms,
            'ltags': self.ltags,
            'lcomment_files': self.lcomment_files,
            'lfiles': self.lfiles,
            'lstruct_numbers': self.lstruct_numbers,
            'lfull_struct_numbers': self.lfull_struct_numbers,
            'lforces': self.lforces,
            'lenergies': self.lenergies
        }

    # OUTCAR to n2p2
    # It's always be a part of input.data
    def write(self,  outfile_name: str=None, append: bool=False):
        """Writes the loaded dataset to an NNP input data file.

        Args:
            outfile_name (str): Output filename, e.g., 'input.data'.
            append (bool, optional): If True, appends to an existing file; otherwise, overwrites it.

        Example:
            >>> data.write("input.data", append=False)
        """
        # remove the file if it exists
        if not append:
            open(outfile_name, 'w').close()

        for atoms, struct_num, full_struct_num, files, tags, comment_files \
                            in zip(self.latoms, self.lstruct_numbers, self.lfull_struct_numbers, self.lfiles, self.ltags, self.lcomment_files):
            self.write_from_ase(outfile_name, atoms, struct_num= struct_num, full_struct_num = full_struct_num, file_name = files, 
                                         tag = tags, comment_file = comment_files)

    # ase to n2p2, only single frame
    def write_from_ase(self, fd: str=None, atoms: Atoms=None, struct_num: int = 1, full_struct_num: int = 1,
                                file_name: str = '', tag: str = 'all', comment_file: str = None, append: bool = True):
        """Writes a single ASE Atoms object to NNP (n2p2) format.

        Args:
            fd (str): Output filename.
            atoms (ase.Atoms): ASE Atoms object containing positions, cell, forces, etc.
            struct_num (int, optional): Current frame index in the source file.
            full_struct_num (int, optional): Total number of frames in the source file.
            file_name (str, optional): Source filename for tracking provenance.
            tag (str, optional): Label tag for the dataset (e.g., 'train', 'test').
            comment_file (str, optional): External comment filename if applicable.
            append (bool, optional): If True, append to existing file.

        Notes:
            This function is typically called by `write()` for batch writing.

        Example:
            >>> data.write_from_ase(fd="input.data", atoms=my_atoms, struct_num=1, full_struct_num=10)
        """
        lattice = np.array(atoms.get_cell())
        # calculator/abc.py get_potential_energy() = calc.get_property('energy', atoms)
        # 'energy' in OUTCAR     : energy(sigma->0) =     -650.294848
        # 'Free energy' in OUTCAR: free  energy   TOTEN  =      -650.29484848 eV
        energy = atoms.get_potential_energy()
        positions = atoms.get_positions()
        # raw forces
        forces = atoms.get_forces(apply_constraint=False)
        chemical_symbols = atoms.get_chemical_symbols()

        # in one outcar file, there are multiple configurations
        # 这里不应设置append, 因为在generate_dataset_from_outcar中已经设置了. 否则当generate中设置'w'时无法输出多帧frames
        # 勘误：上层函数的append参数不要传到这个函数里就ok.
        if not append:
            fd = open(fd, 'w')
        else:
            fd = open(fd, "a")

        # write to fd
        fd.write("begin\n")
        if comment_file not in [None, '']:
            file_name = comment_file
        fd.write("comment | tag={0:s} | frame={1:d}/{2:d} | file={3:s}\n"
                .format(tag, struct_num, full_struct_num, file_name))
        fmt1 = "18.9f"
        # {{}} for keeping {}. {0} is for the first argument - fmt1
        format1_str = "lattice {{0:{0}}} {{1:{0}}} {{2:{0}}}\n".format(fmt1)
        #format1_str = f"lattice {{0:{fmt1}}} {{1:{fmt1}}} {{2:{fmt1}}}\n"
        fmt2 = "12.6f"
        fmt3 = ">5s"
        format2_str = "atom {{0:{0}}} {{1:{0}}} {{2:{0}}} {{3:{1}}} {{4:{1}}} {{5:{1}}} {{6:{0}}} {{7:{0}}} {{8:{0}}}\n".format(fmt2, fmt3)
        fd.write(format1_str.format(lattice[0][0], lattice[0][1], lattice[0][2]))
        fd.write(format1_str.format(lattice[1][0], lattice[1][1], lattice[1][2]))
        fd.write(format1_str.format(lattice[2][0], lattice[2][1], lattice[2][2]))
        for p, f, s in zip(positions, forces, chemical_symbols):
            fd.write(format2_str.format(p[0], p[1], p[2], s, "0.0", "0.0", f[0], f[1], f[2]))
        fd.write(f"energy {{0:{fmt1}}}\n".format(energy))
        fd.write("charge {0:>10s}\n".format("0.0"))
        fd.write("end\n")

    def get_dict(self) -> dict:
        return self.dict

    @staticmethod
    def read_tags_from_datafile(file_path: str = None, default_tag: str = 'all') -> list:
        """Lightweight reader: one ``tag=`` label per structure from input.data.

        Scans only the ``begin`` / ``comment`` / ``end`` lines and the regex
        ``tag=...`` in the comment, without building ASE Atoms objects. About
        30x faster than ``load_from_datafile`` when only the tags are needed
        (e.g. mapping the structure index in trainpoints/trainforces to a tag);
        results are identical, and the order matches the file/structure index.

        Args:
            file_path (str): Path to the input.data file.
            default_tag (str): Tag for structures whose comment has no tag=
                field. Defaults to 'all' (same convention as load_from_datafile).

        Returns:
            list: Tag strings, one per structure, in file order.
        """
        ltag = []
        tag = default_tag
        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                fields = line.split()
                if not fields:
                    continue
                if fields[0] == 'begin':
                    tag = default_tag
                elif fields[0] == 'comment':
                    m = re.search(r'tag=([^\s]+)', line)
                    if m:
                        tag = m.group(1)
                elif fields[0] == 'end':
                    ltag.append(tag)
        return ltag


def _is_orthogonal_cell(cell: np.ndarray, atol: float = 1e-6) -> bool:
    """Return whether the three cell vectors are mutually orthogonal."""
    arr = np.asarray(cell, dtype=float)
    if arr.shape != (3, 3):
        return False
    lengths = np.linalg.norm(arr, axis=1)
    if np.any(lengths <= 0):
        return False
    dots = arr @ arr.T
    off_diag = dots - np.diag(np.diag(dots))
    scale = float(np.max(lengths) ** 2)
    return bool(np.all(np.abs(off_diag) <= atol * scale))


def _stretch_lattice_values_for_epoch_convention(
        lat: str, extr_rvec: np.ndarray, cell_ref: np.ndarray) -> tuple:
    """Convert DFT stretch row-vector lengths to the LAMMPS epoch-scan convention.

    The VASP archive stores FCC/BCC equilibrium structures as primitive cells, so
    ``my_read_stretch`` returns primitive row-vector lengths. LAMMPS y_post scans
    use conventional cubic cells for FCC/BCC. HCP already uses the same hexagonal
    ``a``/``c`` convention on both sides.
    """
    rvec = np.asarray(extr_rvec, dtype=float)
    a = float(rvec[0])
    c = float(rvec[2])
    if _is_orthogonal_cell(cell_ref):
        return a, c
    if lat == 'fcc':
        factor = np.sqrt(2.0)
        return a * factor, c * factor
    if lat == 'bcc':
        factor = 2.0 / np.sqrt(3.0)
        return a * factor, c * factor
    return a, c


# DFT (VASP) 参考基线，作为 PeiN2p2.post_epoch_scan 的对照量：post_epoch_scan 逐 epoch 读取
# LAMMPS (pair_style hdnnp) 的 p_post_stretch.txt / y_post_cij_energy.txt / y_post_gsfe.txt，
# 本函数用同一批读取器（my_read_stretch / read_cij_energy / read_output）从 construct_dataset 的
# DFT 计算归档里读出同样的物理量，键名与 post_epoch_scan 的行字典一一对应，便于画图时叠加 DFT 水平参考线。
# 相<->标签映射（construct_dataset/README.md）：A11=HCP, A12=FCC, A13=BCC；A21=HCP 缺陷, A22=FCC 缺陷。
#   stretch  {tag}/y_stretch/p_post_stretch.txt        -> a_<lat>, c_<lat>, E_<lat>, ca_hcp, dE_*_fcc
#              FCC/BCC 的 VASP primitive rvector 会换算为 LAMMPS y_post 的 conventional cubic a
#   cij      {tag}/{cij_subdir}/y_post_cij_energy.txt   -> C11_<lat> ... C44_<lat>
#              默认 cij_subdir='y_cij_energy_small'（小应变，harmonic 弹性常数），与训练/后处理统一
#   gsfe     {tag}/{type}/y_gsfe_{type}/y_post_gsfe.txt -> usf_<type>=max(gamma), sf_<type>=最后一个 gamma
# 缺文件不报错，跳过该项（DFT 归档可能只覆盖部分相/滑移系），返回的子字典只含读到的键。
def read_dft_reference(dir_dft_root,
                       lats=('fcc', 'bcc', 'hcp'),
                       cij_keys=('C11', 'C12', 'C13', 'C33', 'C44'),
                       stretch_tag=None, cij_tag=None, cij_subdir='y_cij_energy_small',
                       gsfe_tag=None, gsfe_types=None,
                       verbose=True) -> dict:
    """Batch-read DFT (VASP) reference properties into one dict.

    Reads the same per-phase stretch / Cij / GSFE post-processing files that
    ``PeiN2p2.post_epoch_scan`` reads per training epoch, but from the DFT
    archive under ``dir_dft_root`` (e.g. ``construct_dataset/calculate``). Uses
    the identical readers (:func:`mymetal.post.stretch.my_read_stretch`,
    :func:`mymetal.post.Cij_energy.read_cij_energy`,
    :func:`mymetal.post.gsfe.read_output`) so the returned keys line up with the
    epoch-scan row dicts and can be overlaid as horizontal reference lines.
    FCC/BCC DFT primitive-cell stretch vectors are converted to conventional
    cubic lattice constants to match LAMMPS y_post.

    Args:
        dir_dft_root (str | Path): Root of the DFT archive holding the ``A1x``/
            ``A2x`` tag directories.
        lats (tuple): Phases to read for stretch/cij. Default ``('fcc','bcc','hcp')``.
        cij_keys (tuple): Elastic-constant keys to pull from each cij file.
        stretch_tag (dict, optional): phase -> tag for stretch files.
            Default ``{'fcc':'A12-1','bcc':'A13-1','hcp':'A11-1'}``.
        cij_tag (dict, optional): phase -> tag for cij files.
            Default ``{'fcc':'A12-2','bcc':'A13-2','hcp':'A11-2'}``.
        cij_subdir (str): Sub-directory under each cij tag holding the
            post-processed elastic-constant file. Default
            ``'y_cij_energy_small'`` (small-strain / harmonic Cij, matching the
            small-strain cij training set and the epoch-scan LAMMPS evaluation).
            Pass ``'y_cij_energy'`` for the legacy large-strain fit.
        gsfe_tag (dict, optional): phase -> tag for gsfe files.
            Default ``{'fcc':'A22-2','hcp':'A21-2'}``.
        gsfe_types (dict, optional): phase -> list of slip-system type names.
            Default matches ``post_epoch_scan``'s ``gsfe_types``.
        verbose (bool): Print a line for each missing/skipped file.

    Returns:
        dict: ``{'stretch': {...}, 'cij': {...}, 'gsfe': {...}}`` with flat
        sub-dicts keyed like the ``post_epoch_scan`` columns (no ``epoch``).
        ``stretch``: ``a_<lat>``, ``c_<lat>``, ``E_<lat>`` (eV/atom) plus derived
        ``ca_hcp`` (-) and ``dE_hcp_fcc`` / ``dE_bcc_fcc`` (meV/atom). FCC/BCC
        ``a``/``c`` are conventional cubic lengths when the DFT archive uses
        primitive VASP cells;
        ``cij``: ``<C..>_<lat>`` (GPa); ``gsfe``: ``usf_<type>`` / ``sf_<type>``
        (mJ/m^2; forced from the GSFE gamma table as max(gamma) / last gamma).
    """
    # 与 post_epoch_scan 一致：就地局部 import，避免给 dataset 顶层引入 mymetal.post 依赖
    from contextlib import redirect_stdout
    from io import StringIO
    from mymetal.post.stretch import my_read_stretch
    from mymetal.post.Cij_energy import read_cij_energy
    from mymetal.post.gsfe import read_output

    if stretch_tag is None:
        stretch_tag = {'fcc': 'A12-1', 'bcc': 'A13-1', 'hcp': 'A11-1'}
    if cij_tag is None:
        cij_tag = {'fcc': 'A12-2', 'bcc': 'A13-2', 'hcp': 'A11-2'}
    if gsfe_tag is None:
        gsfe_tag = {'fcc': 'A22-2', 'hcp': 'A21-2'}
    if gsfe_types is None:
        gsfe_types = {'fcc': ['FCC_100', 'FCC_111'],
                      'hcp': ['HCP_basal', 'HCP_prism1w', 'HCP_pyr1w', 'HCP_pyr2']}

    root = Path(dir_dft_root)

    # 1) stretch：各相平衡 a/c (Extr rvector) 与 E/atom (Extr infos)
    stretch = {}
    for lat in lats:
        tag = stretch_tag.get(lat)
        if tag is None:
            continue
        ps = root / tag / 'y_stretch' / 'p_post_stretch.txt'
        if not ps.is_file():
            if verbose:
                print(f'skip dft stretch {lat}: missing {ps}')
            continue
        with redirect_stdout(StringIO()):              # my_read_stretch 会回显注释头，抑制
            *_, cell_ref, (_, _, extr_y, extr_rvec) = my_read_stretch(str(ps))
        a_lat, c_lat = _stretch_lattice_values_for_epoch_convention(lat, extr_rvec, cell_ref)
        stretch[f'a_{lat}'] = a_lat
        stretch[f'c_{lat}'] = c_lat
        stretch[f'E_{lat}'] = float(extr_y)
    if 'a_hcp' in stretch and 'c_hcp' in stretch:
        stretch['ca_hcp'] = stretch['c_hcp'] / stretch['a_hcp']
    if 'E_hcp' in stretch and 'E_fcc' in stretch:
        stretch['dE_hcp_fcc'] = (stretch['E_hcp'] - stretch['E_fcc']) * 1e3
    if 'E_bcc' in stretch and 'E_fcc' in stretch:
        stretch['dE_bcc_fcc'] = (stretch['E_bcc'] - stretch['E_fcc']) * 1e3

    # 2) cij：各相 C11..C44（立方相文件同样写满 5 个，C13/C33 与 C12/C11 近似重复）
    cij = {}
    for lat in lats:
        tag = cij_tag.get(lat)
        if tag is None:
            continue
        pc = root / tag / cij_subdir / 'y_post_cij_energy.txt'
        if not pc.is_file():
            if verbose:
                print(f'skip dft cij {lat}: missing {pc}')
            continue
        d = read_cij_energy(str(pc))
        for k in cij_keys:
            cij[f'{k}_{lat}'] = d[k]

    # 3) gsfe：逐滑移系读取 usf=max(gamma)、sf=最后一个 gamma
    gsfe = {}
    for phase, ltype in gsfe_types.items():
        tag = gsfe_tag.get(phase)
        if tag is None:
            continue
        for t in ltype:
            pg = root / tag / t / f'y_gsfe_{t}' / 'y_post_gsfe.txt'
            if not pg.is_file():
                if verbose:
                    print(f'skip dft gsfe {t}: missing {pg}')
                continue
            g = read_output(str(pg))
            gsfe[f'usf_{t}'] = g['usf_max'] if g['usf_max'] is not None else np.nan
            gsfe[f'sf_{t}'] = g['sf_min'] if g['sf_min'] is not None else np.nan

    return {'stretch': stretch, 'cij': cij, 'gsfe': gsfe}
