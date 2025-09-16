"""
dataset module

This module provides functions for creating, reading, and writing datasets in NNP format.
It includes functions to extract atomic configurations, energies, and forces from VASP OUTCAR files,
as well as to write these configurations into a dataset file compatible with NNP tools.

Functions:
    - write_nnp_data_from_ase: Writes atomic configuration data from an ASE Atoms object to a file in NNP format.
    - generate_dataset_from_outcar: Extracts atomic structure frames from a VASP OUTCAR file and writes them to a dataset file.
    - generate_dataset_from_outcar_n2p2: Parses a VASP OUTCAR file and converts one or more structures into NNP input format.
    - read_nnp_dataset: Parses a NNP-format dataset file into a structured dictionary of ASE Atoms objects grouped by tag.
    - generate_dataset_from_dict: Writes ASE Atoms structures from a grouped dictionary into a NNP-compatible dataset file.
"""


from ase.calculators.singlepoint import SinglePointCalculator
from ase.io.vasp import read_vasp_out
from ase import Atoms
import os
import numpy as np
import pandas as pd
import re
from collections import defaultdict

# ase to n2p2, only single frame
def write_nnp_data_from_ase(fd: str=None, atoms: Atoms=None, index: int = 1, largest_index: int = 1,
                            file_name: str = '', tag: str = 'all', append: bool = True):
    """
    Writes atomic configuration data from an ASE Atoms object to a file in NNP format.

    Args:
        fd (str): Path to the output file to append data.
        atoms (Atoms): ASE Atoms object containing atomic structure, positions, and forces.
        index (int): Index of the current structure frame (1-based).
        largest_index (int): Total number of frames in the dataset.
        file_name (str): Source OUTCAR file name, used in the comment line.
        tag (str): Custom tag for labeling the structure origin or category.
        append (bool): If True, appends to the file; if False, overwrites it.
        
    Writes:
        A single structure block including lattice vectors, atom positions, forces,
        energy, and charge to the specified file.

    Note:
        The force, energy is extracted from the calculator object. 
        For force, the apply_constraint is set to False to get the raw forces.
        Foe energy, the potential energy is energy(sigma->0).
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
    if append:
        fd = open(fd, "a")
    else:
        fd = open(fd, "w")
    # write to fd
    # "fd" is actually a file-descriptor thanks to @writer
    fd.write("begin\n")
    fd.write("comment | tag={0:s} | structure_number={1:d}/{2:d} frames | source_file_name={3:s}\n"
             .format(tag, index, largest_index, file_name))
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

# OUTCAR to n2p2
def generate_dataset_from_outcar(outcarfile: str=None, outfile_name: str=None, append: bool=False,
                                    index: str = ':', tag: str = 'all'):
    """
    Extracts one or more atomic structure frames from a VASP OUTCAR file and writes them to a dataset file.

    Args:
        outcarfile (str): Path to the OUTCAR file to read atomic structures from.
        outfile_name (str): Path to the output file where data will be written.
        append (bool): If False, clears the output file before writing; if True, appends to it.
        index (str): Frame selection string, e.g., ':' for all frames or '-1' for the last.
        tag (str): Label to identify the source or purpose of the data (used in the output comment).

    Writes:
        One or more frames of atomic data to the output file in NNP-compatible format.

    Example:
        >>> file_name = './workdir2/OUTCAR'
        >>> outfile_name = './input-test.data'
        >>> my_generate_dataset_from_outcar(file_name, outfile_name, tag='elastic')
    """
    # Read whole atoms configurations from OUTCAR
    # index=':' means read all configurations
    outcar = read_vasp_out(outcarfile, index=index)
    num = len(outcar)
    # remove the file if it exists
    if not append == True:
        open(outfile_name, 'w').close()
    for i, atom in enumerate(outcar):
        # use 'a' to append all configurations
        # in one outcar file, there are multiple configurations
        write_nnp_data_from_ase(outfile_name, atom, i+1, num, outcarfile, tag)

# taken from n2p2 tools
# The code structure remains unchanged, only encapsulated as a function
def generate_dataset_from_outcar_n2p2(file_name: str=None, outfile_name: str=None):
    """
    Parses a VASP OUTCAR file and converts one or more structures into NNP input format.

    This function reads atomic configurations, lattice vectors, total energies, and atomic
    forces from a VASP OUTCAR file, including cases with multiple ionic steps (e.g., MD or relaxation).
    The extracted data is written to a file in the NNP-compatible "input.data" format.

    Args:
        file_name (str): Path to the OUTCAR file generated by VASP.
        outfile_name (str): Path to the output file where parsed data will be written. 
                            If None, output is printed to stdout.

    Raises:
        RuntimeError: If the number of extracted lattice blocks, energies, and atomic data blocks
                      do not match (indicating a malformed or incomplete OUTCAR file).

    Example:
        >>> file_name = './workdir2/OUTCAR'
        >>> outfile_name = './input2.data'
        >>> generate_dataset_from_outcar(file_name, outfile_name)
    """
    # Read in the whole file first.
    f = open(file_name, "r")
    lines = [line for line in f]
    f.close()

    # If OUTCAR contains ionic movement run (e.g. from an MD simulation) multiple
    # configurations may be present. Thus, need to prepare empty lists.
    lattices   = []
    energies   = []
    atom_lists = []

    # Loop over all lines.
    elements = []
    for i in range(len(lines)):
        line = lines[i]
        # Collect element type information, expecting VRHFIN lines like this:
        #
        # VRHFIN =Cu: d10 p1
        #
        if "VRHFIN" in line:
            elements.append(line.split()[1].replace("=", "").replace(":", ""))
        # VASP specifies how many atoms of each element are present, e.g.
        #
        # ions per type =              48  96
        #
        if "ions per type" in line:
            atoms_per_element = [int(it) for it in line.split()[4:]]
        # Simulation box may be specified multiple times, I guess this line
        # introduces the final lattice vectors.
        # Simulation box may be specified multiple times, I guess this line
        # introduces the final lattice vectors.
        # Here has bug: 
        #       direct lattice vectors                 reciprocal lattice vectors
        #    -3.660342136 -3.660342136  0.000000000    -0.136599253 -0.136599253  0.000000000
        #    10.981026407-10.981026407  0.000000000     0.045533084 -0.045533084  0.000000000
        #     0.000000000  0.000000000 41.982653687     0.000000000  0.000000000  0.023819361
        #    10.981026407-10.981026407 is not correct, should be -10.981026407 -10.981026407
        if "VOLUME and BASIS-vectors are now" in line:
            lattices.append([lines[i+j].split()[0:3] for j in range(5, 8)])
        # Total energy is found in the line with "energy  without" (2 spaces) in
        # the column with sigma->0:
        #
        # energy  without entropy=     -526.738461  energy(sigma->0) =     -526.738365
        #
        if "energy  without entropy" in line:
            energies.append(line.split()[6])
        # Atomic coordinates and forces are found in the lines following
        # "POSITION" and "TOTAL-FORCE".
        if "POSITION" in line and "TOTAL-FORCE" in line:
            atom_lists.append([])
            count = 0
            for ei in range(len(atoms_per_element)):
                for j in range(atoms_per_element[ei]):
                    atom_line = lines[i+2+count]
                    atom_lists[-1].append(atom_line.split()[0:6])
                    atom_lists[-1][-1].extend([elements[ei]])
                    count += 1

    # Sanity check: do all lists have the same length. 
    if not (len(lattices) == len(energies) and len(energies) == len(atom_lists)):
        raise RuntimeError("ERROR: Inconsistent OUTCAR file.")

    # Open output file or write to stdout.
    if outfile_name is not None:
        f = open(outfile_name, "w")
    else:
        f = sys.stdout

    # Write configurations in "input.data" format.
    for i, (lattice, energy, atoms) in enumerate(zip(lattices, energies, atom_lists)):
        f.write("begin\n")
        f.write("comment source_file_name={0:s} structure_number={1:d}\n".format(file_name, i + 1))
        f.write("lattice {0:s} {1:s} {2:s}\n".format(lattice[0][0], lattice[0][1], lattice[0][2]))
        #print(type(lattice[0][0]))
        f.write("lattice {0:s} {1:s} {2:s}\n".format(lattice[1][0], lattice[1][1], lattice[1][2]))
        f.write("lattice {0:s} {1:s} {2:s}\n".format(lattice[2][0], lattice[2][1], lattice[2][2]))
        for a in atoms:
            f.write("atom {0:s} {1:s} {2:s} {3:s} {4:s} {5:s} {6:s} {7:s} {8:s}\n".format(a[0], a[1], a[2], a[6], "0.0", "0.0", a[3], a[4], a[5]))
        f.write("energy {0:s}\n".format(energy))
        f.write("charge {0:s}\n".format("0.0"))
        f.write("end\n")

# n2p2 to dict (include the atoms object)
def read_nnp_dataset(file_path: str= None) -> dict:
    """
    Parses a NNP-format dataset file into a structured dictionary of ASE Atoms objects grouped by tag.

    This function reads a dataset file (e.g. 'input.data') formatted using NNP tools. It extracts atomic 
    structures, energy, force, and metadata such as tag, source file, and structure number.

    Args:
        file_path (str): Path to the input dataset file in NNP format.

    Returns:
        dict: A dictionary of structure groups, where each tag maps to a dictionary containing:
            - 'atoms' (List[ase.Atoms]): List of ASE Atoms objects.
            - 'source_file' (List[str]): Source OUTCAR file names.
            - 'struct_number' (List[int]): Structure index within the source file.
            - 'full_struct_number' (List[int]): Total number of structures in the source file.

            Example structure:
            {
                "train": {
                    "atoms": [Atoms1, Atoms2, ...],
                    "source_file": ["OUTCAR", "OUTCAR", ...],
                    "struct_number": [1, 2, ...],
                    "full_struct_number": [100, 100, ...]
                },
                ...
            }

    Raises:
        ValueError: If the number of positions, forces, and symbols does not match.
        ValueError: If the input file path is None.

    Example:
        >>> dataset = read_nnp_dataset("input.data")
        >>> len(dataset["train"]["atoms"])
        10
    """
    if file_path is None:
        raise ValueError("The input file path is None.")
    
    with open(file_path, 'r') as f:
        lines = f.readlines()
    i = 0
    n = len(lines)

    dict = defaultdict(lambda: {'atoms': [], 'source_file': [], 
                                'struct_number': [], 'full_struct_number': []})

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
                    # Example format: comment | tag=xx | structure_number=1/12 frames | source_file_name=xxx
                    comment = line
                    tag_match = re.search(r'tag=([^\s]+)', comment)
                    struct_match = re.search(r'structure_number=(\d+)(?:/(\d+))?', comment)
                    file_match = re.search(r'source_file_name=([^\s]+)', comment)
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
            dict[tag]['struct_number'].append(struct_number)
            dict[tag]['full_struct_number'].append(full_struct_number)
            dict[tag]['source_file'].append(source_file)
            dict[tag]['atoms'].append(atoms)
        i = i+1
    return dict

# dict to n2p2
def generate_dataset_from_dict(dict: dict={}, outfile_name: str = 'input.data', append: bool = False):
    """
    Writes ASE Atoms structures from a grouped dictionary into a NNP-compatible dataset file.

    This function takes the dictionary output of `read_nnp_dataset()` and writes its contents
    into a dataset file in NNP input format using `write_nnp_data_from_ase`.

    Args:
        dict (dict): A dictionary grouped by tag, containing lists of atoms and metadata.
                     Should follow the structure returned by `read_nnp_dataset()`.
        outfile_name (str): Path to the output file. Defaults to 'input.data'.
        append (bool): If True, appends to the file; otherwise, overwrites. Defaults to False.

    Raises:
        ValueError: If the input dictionary is None.

    Example:
        >>> dict = read_nnp_dataset("input.data")
        >>> generate_dataset_from_dict(dict, "./output.data", append=False)
    """
    if dict is None:
        raise ValueError("The input dictionary is None.")
    if not append == True:
        open(outfile_name, 'w').close()
    for key in dict.keys():
        keydict = dict[key]
        atoms_list = keydict['atoms']
        source_file_list = keydict['source_file']
        struct_number_list = keydict['struct_number']
        full_struct_number_list = keydict['full_struct_number']
        for a, sf, sn, fsn in zip(atoms_list, source_file_list, struct_number_list, full_struct_number_list):
            write_nnp_data_from_ase(outfile_name, a, sn, fsn, sf, key)