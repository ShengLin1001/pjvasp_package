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

Change Log:
    - 2025.10.19 Integrated all functions into the nnpdata class for better encapsulation.
"""


from ase.calculators.singlepoint import SinglePointCalculator
from ase.io.vasp import read_vasp_out
from ase import Atoms
import os
import numpy as np
import pandas as pd
import re
from collections import defaultdict


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

        self.outcar_dict = {
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