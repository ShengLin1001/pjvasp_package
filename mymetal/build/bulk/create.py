"""
create submodule

This submodule provides functions for creating various (hkl)-oriented crystal cells,
replacing the previous approach of using `generate_bulk_from_film()` in the film module.

Functions:
    - create_fcc_111: Creates an FCC (111) plane structure.
    - create_hcp_basal: Creates an HCP (0001) basal plane structure.
    - create_hcp_prism1: Creates an HCP (10-10) prism
    - create_hcp_prism2: Creates an HCP (10-11) prism.

Background:
    Crystal structures and their characteristic dislocations and slip systems are summarized below.

    Full dislocations:
        - fcc: a/2 <110>, b = 1/2 <110>
        - bcc: a/2 <111>, b = 1/2 <111>
        - hcp: a/3 <11-20>, b = 1/3 <11-20>

    Major slip planes:
        - fcc: (111)
        - bcc: {110}, {112}, {123}
        - hcp: (0001), (10-10), (10-11)

    Partial dislocations:
        - fcc: a/6 <112> on (111) — Shockley partial
        - bcc: a/2[11-1] → a/3[11-1] + a/6[11-1]
        - hcp: a/3[11-20] → 1/3[10-10] + 1/3[01-10]

"""

from ase.lattice.cubic import FaceCenteredCubic, BodyCenteredCubic
from ase.lattice.hexagonal import HexagonalClosedPacked
from mymetal.universal.atom.moveatom import move_atoms
from ase import Atoms
import numpy as np
from myvasp import vasp_func as vf

INPLANE_SHIFT = 1e-5  # to avoid atoms sitting exactly on the xlo/ylo face of the (tilted) cell

# FCC (111) Plane
def create_fcc_111(a: float = None, size: tuple = (1, 1, 1),
                   pbc: tuple = (1, 1, 1), symbol: str = 'Au') -> Atoms:
    direction = [[-1, 1, 0],
                    [0, -1, 1],
                    [1, 1, 1]]
    miller = [None, None, [1, 1, 1]]
    atoms = FaceCenteredCubic(directions=direction,
                            miller = miller,
                            size=size, symbol=symbol, pbc=pbc,
                            latticeconstant = a)
    # the interlayer spacing of fcc (111) is a/sqrt(3) ~ 1.732 a
    atoms = move_atoms(atoms, [INPLANE_SHIFT, INPLANE_SHIFT, 0.1*a], if_scale_position=False)
    atoms.wrap()
    return atoms

# HCP (0001) Basal Plane
def create_hcp_basal(a: float = None, c: float = None, size: tuple = (1, 1, 1),
                   pbc: tuple = (1, 1, 1), symbol: str = 'Au') -> Atoms:
    direction = [[2, -1, -1, 0], 
        [-1, 2, -1, 0],
        [ 0, 0, 0,  1]]
    miller = [None, None, None] # hcp can't set miller index
    atoms = HexagonalClosedPacked(directions=direction,
                                         miller = miller,
                                         size=size, symbol=symbol, pbc=pbc,
                                         latticeconstant = {'a': a, 'c': c})
    # the interlayer spacing of hcp (0001) is c/2
    atoms = move_atoms(atoms, [0,0,0.1*c], if_scale_position=False)
    # nudge in-plane so no atom sits exactly on the xlo/ylo face of the (tilted) cell.
    # otherwise the corner atom at (0,0) gets pushed just outside the box by float noise
    # during a LAMMPS minimize, and write/read_restart drops it
    # ("Did not assign all restart atoms correctly"). Pure translation -> physics unchanged.
    atoms = move_atoms(atoms, [INPLANE_SHIFT, INPLANE_SHIFT, 0.0], if_scale_position=True)
    atoms.wrap()
    return atoms

# HCP (10-10) Prism I, Wide or Narrow
def create_hcp_prism1(a: float = None, c: float = None, size: tuple = (1, 1, 1),
                   pbc: tuple = (1, 1, 1), symbol: str = 'Au',
                    mode: str = 'wide') -> Atoms:
    direction = [[-1, 2, -1, 0],
                [ 0, 0, 0, 1],
                [1, 0, -1, 0]]
    miller = [None, None, None]
    atoms = HexagonalClosedPacked(directions=direction,
                                            miller = miller,
                                            size=size, symbol=symbol, pbc=pbc,
                                            latticeconstant = {'a': a, 'c': c})
    # Now, the sequence: narrow, wide, narrow, wide, ...
    # Interplanar spacing of wide: 1/sqrt(3) * a ~ 0.577 * a
    # Interplanar spacing of narrow: 1/(2*sqrt(3)) * a ~ 0.289 * a
    if mode == 'wide':
        atoms = move_atoms(atoms, [INPLANE_SHIFT, INPLANE_SHIFT, 0.1*a], if_scale_position=False)
    elif mode == 'narrow':
        atoms = move_atoms(atoms, [INPLANE_SHIFT, INPLANE_SHIFT, -0.05*a], if_scale_position=False)
    atoms.wrap()
    return atoms

# HCP (10-11) Prism II
def create_hcp_prism2(a: float = None, c: float = None, size: tuple = (1, 1, 1),
                   pbc: tuple = (1, 1, 1), symbol: str = 'Au') -> Atoms:
    direction = [[-1, 1, 0, 0],
                [ 0, 0, 0, 1],
                [1 , 1, -2, 0]]
    miller = [None, None, None]
    atoms = HexagonalClosedPacked(directions=direction,
                                            miller = miller,
                                            size=size, symbol=symbol, pbc=pbc,
                                            latticeconstant = {'a': a, 'c': c})
    # the interlayer spacing of hcp (10-11) is 0.5 * a
    atoms = move_atoms(atoms, [INPLANE_SHIFT, INPLANE_SHIFT, 0.1*a], if_scale_position=False)
    atoms.wrap()
    return atoms

# todo: BCC (110) Plane


# Vectorized drop-in replacement for myvasp.vasp_create.create_supercell.
# The original builds atoms_pos by np.vstack inside a 4-deep Python loop, which
# is O(N^2) in memory copies and dominates the runtime of vasp_create_hcp_pyr1/
# pyr2 (and the other create_* helpers) once the cell holds thousands of atoms.
# This version preallocates and broadcasts -> O(N), while keeping the exact same
# atom ordering (cell index k-slow,j,i-fast outermost; motif index m innermost)
# so the produced model is identical to vf.create_supercell.
def create_supercell_fast(latt: np.ndarray, motif: np.ndarray,
                          ncell) -> Atoms:
    """Build a supercell by tiling a primitive lattice + motif basis.

    Vectorized drop-in replacement for ``myvasp.vasp_create.create_supercell``.
    Preallocates and broadcasts instead of vstack inside a 4-deep Python loop,
    reducing runtime from O(N²) memory copies to O(N) while keeping the exact
    same atom ordering (cell index k-slow, j, i-fast outermost; motif index m
    innermost) as the original.

    Args:
        latt (np.ndarray): 3×3 lattice matrix, row vectors (Å).
        motif (np.ndarray): (nmotif, 3) fractional basis positions.
        ncell (array-like): (nx, ny, nz) tiling factors.

    Returns:
        ase.Atoms: supercell with PBC=[1,1,1] and chemical symbols set to
        the default (element 1). The caller is responsible for assigning
        real symbols.
    """
    latt = np.asarray(latt, dtype=float)
    motif = np.asarray(motif, dtype=float)
    nx, ny, nz = (int(ncell[0]), int(ncell[1]), int(ncell[2]))

    # cell-origin indices in the original loop order: k slowest, then j, then i
    K, J, I = np.meshgrid(np.arange(nz), np.arange(ny), np.arange(nx),
                          indexing='ij')
    I, J, K = I.ravel()[:, None], J.ravel()[:, None], K.ravel()[:, None]
    # accumulate elementwise in the same ((i*l0 + j*l1) + k*l2) order as the
    # original loop, so there is no floating-point re-association difference
    refp = I * latt[0, :] + J * latt[1, :] + K * latt[2, :]     # (ncells, 3)

    motif_cart = motif @ latt                                  # (nmotif, 3)

    # broadcast: cell index slow, motif index fast -> matches original ordering
    atoms_pos = (refp[:, None, :] + motif_cart[None, :, :]).reshape(-1, 3)

    superlatt = latt.copy()
    for i in np.arange(3):
        superlatt[i, :] = superlatt[i, :] * ncell[i]

    atoms = Atoms(cell=superlatt, positions=atoms_pos, pbc=[1, 1, 1])
    natoms = atoms.positions.shape[0]
    atoms.set_chemical_symbols(np.ones([natoms, 1]))
    return atoms


##### This part is taken from
# https://github.com/BinglunYin/myalloy_package/blob/master/myvasp/vasp_create_hcp.py

#===================
# examples:

# a = 3.23415
# ca = 1.5992

## For pry1 and prism1:
##      bp: 33 for W, -33 for N

# vasp_create_hcp_basal(a, ca, np.array([1, 1, 10]) )
# vasp_create_hcp_prism1(a, ca, np.array([1, 1, 16]), bp = 33)
# vasp_create_hcp_pyr1(a, ca, np.array([1, 1, 16]), bp=-33)
# vasp_create_hcp_pyr1(a, ca, np.array([1, 1, 16]), bp=33)
# vasp_create_hcp_pyr2(a, ca, np.array([1, 1, 10]) )

def vasp_create_hcp_basal(a, ca, ncell, bp=33):
    """Build an HCP (0001) basal-plane supercell via lattice + motif.

    Uses the primitive hexagonal lattice (a, sqrt(3)/2, ca) with a 2-atom
    motif at (0,0,0) and (1/3, 2/3, 1/2). The ``bp`` flag shifts atoms
    in +z by 0.1*a (bp=33) or not (bp=0) to avoid surface atoms sitting
    on the boundary.

    Args:
        a (float): HCP lattice constant a (Å).
        ca (float): c/a ratio (dimensionless); c = ca * a.
        ncell (np.ndarray): (nx, ny, nz) tiling factors.
        bp (int): boundary-position flag; 33 shifts +z, 0 does not.

    Returns:
        ase.Atoms: HCP basal supercell. ``atoms.pos_a0`` stores the
        reference lattice constant.

    Note:
        Depends on the optional ``myvasp`` companion package.
    """
    print('==> create hcp basal plane: ')
    print(a, ca, ncell, bp)

    latt = np.array([
        [1.0, 0, 0],
        [-0.5, np.sqrt(3)/2, 0],   
        [0, 0, ca],
    ]) * a
    
    motif = np.array([
        [0, 0, 0],
        [1/3,  2/3,  1/2],
    ])

    atoms = create_supercell_fast(latt, motif, ncell)

    if bp == 33:
        atoms.positions = atoms.positions + np.array([0, 0, 0.1])*a
        atoms.wrap()

    atoms.pos_a0 = a 
    return atoms
    #vf.my_write_vasp(atoms, filename='POSCAR', vasp5=True)

def vasp_create_hcp_basal_ortho(a, ca, ncell, bp=33):
    """Build an orthorhombic HCP (0001) basal-plane supercell.

    Uses an orthorhombic lattice (a, sqrt(3)*a, ca*a) with a 4-atom motif,
    equivalent to two primitive cells. Suitable for LAMMPS which prefers
    orthogonal boxes.

    Args:
        a (float): HCP lattice constant a (Å).
        ca (float): c/a ratio (dimensionless).
        ncell (np.ndarray): (nx, ny, nz) tiling factors.
        bp (int): boundary-position flag; 33 shifts +z, 0 does not.

    Returns:
        ase.Atoms: orthorhombic HCP basal supercell.

    Note:
        Depends on the optional ``myvasp`` companion package.
    """
    print('==> create hcp basal_ortho plane: ')
    print(a, ca, ncell, bp)

    latt = np.array([
        [ 1.0,          0,  0],
        [   0, np.sqrt(3),  0],   
        [   0,          0, ca],
    ]) * a
    
    motif = np.array([
        [0.0,    0,    0],
        [1/2,  1/2,    0],
        [  0,  2/6,  1/2],
        [1/2,  5/6,  1/2],
    ])

    atoms = create_supercell_fast(latt, motif, ncell)
    
    if bp == 33:
        atoms.positions = atoms.positions + np.array([0, 0, 0.1])*a
        atoms.wrap()

    atoms.pos_a0 = a 
    return atoms
    #vf.my_write_vasp(atoms, filename='POSCAR', vasp5=True)

def vasp_create_hcp_prism1(a, ca, ncell, bp=33):
    """Build an HCP (10-10) prism I supercell.

    Starts from the basal cell, then applies ``make_SFP_xy`` (surface
    frame projection) and ``make_a3_ortho`` to reorient the cell so
    that a3 points along the prism direction. The ``bp`` flag selects
    wide (bp=33, +0.1*a shift) or narrow (bp=-33, -0.05*a shift).

    Args:
        a (float): HCP lattice constant a (Å).
        ca (float): c/a ratio (dimensionless).
        ncell (np.ndarray): (nx, ny, nz) tiling factors for the *basal*
            cell; internally swapped for prism orientation.
        bp (int): 33 = wide, -33 = narrow, 0 = no shift.

    Returns:
        ase.Atoms: HCP prism I supercell.

    Note:
        Depends on the optional ``myvasp`` companion package.
    """
    print('==> create hcp prism plane: ')
    print(a, ca, ncell, bp)

    ncell2 = ncell.copy()
    for i in np.arange(3):
        ncell2[i] = ncell[ np.mod(i+2, 3) ]
    
    atoms = vasp_create_hcp_basal(a, ca, ncell2, bp=0)

    #atoms = vf.my_read_vasp('POSCAR')
    atoms = vf.make_SFP_xy(atoms, i1=1)
    atoms = vf.make_a3_ortho(atoms)

    if bp == 33:
        print('==> create prism1-W')
        atoms.positions = atoms.positions + np.array([0, 0, 0.1])*a
        atoms.wrap()
    
    elif bp == -33:
        print('==> create prism1-N')
        atoms.positions = atoms.positions + np.array([0, 0, -0.1])*a
        atoms.wrap()
    
    atoms.pos_a0 = a 
    return atoms
    #vf.my_write_vasp(atoms, filename='POSCAR', vasp5=True)

def vasp_create_hcp_pyr1(a, ca, ncell, bp=33):
    """Build an HCP (10-11) pyramidal I supercell.

    Uses a triclinic lattice with a3 tilted to mix basal and c-directions.
    Applies ``make_SFP_xy`` and ``make_a3_ortho`` after construction.
    The ``bp`` flag selects wide (bp=33, +0.1*a) or narrow (bp=-33, -0.05*a).

    Args:
        a (float): HCP lattice constant a (Å).
        ca (float): c/a ratio (dimensionless).
        ncell (np.ndarray): (nx, ny, nz) tiling factors.
        bp (int): 33 = wide, -33 = narrow, 0 = no shift.

    Returns:
        ase.Atoms: HCP pyramidal I supercell.

    Note:
        Depends on the optional ``myvasp`` companion package.
    """
    print('==> create hcp basal pyr1: ')
    print(a, ca, ncell, bp)

    latt = np.array([
        [ 1.0, 0, 0],
        [ 0.5, np.sqrt(3)/2, 0],   
        [-1.0, 0, ca],
    ]) * a
    
    motif = np.array([
        [0, 0, 0],
        [1/6, 2/3, 1/2],
    ])

    ncell2 = ncell.copy()
    for i in np.arange(3):
        ncell2[i] = ncell[ np.mod(i+2, 3) ]
   
    atoms = create_supercell_fast(latt, motif, ncell2)
    atoms = vf.make_SFP_xy(atoms, i1=1)
    atoms = vf.make_a3_ortho(atoms)

    if bp == 33:
        print('==> create pry1-W')
        atoms.positions = atoms.positions + np.array([0, 0, 0.1])*a
        atoms.wrap()

    elif bp == -33:
        print('==> create pry1-N')
        atoms.positions = atoms.positions + np.array([0, 0, -0.05])*a
        atoms.wrap()
    
        
    atoms.pos_a0 = a 
    return atoms
    # vf.my_write_vasp(atoms, filename='POSCAR', vasp5=True)

def vasp_create_hcp_pyr2(a, ca, ncell, bp=33):
    """Build an HCP (10-12) pyramidal II supercell.

    Uses an orthorhombic base lattice with a3 tilted. 4-atom motif.
    Applies ``make_SFP_xy`` and ``make_a3_ortho`` after construction.

    Args:
        a (float): HCP lattice constant a (Å).
        ca (float): c/a ratio (dimensionless).
        ncell (np.ndarray): (nx, ny, nz) tiling factors.
        bp (int): 33 shifts +z, 0 does not.

    Returns:
        ase.Atoms: HCP pyramidal II supercell.

    Note:
        Depends on the optional ``myvasp`` companion package.
    """
    print('==> create hcp basal pyr2: ')
    print(a, ca, ncell, bp)

    latt = np.array([
        [ 1.0, 0, 0],
        [ 0.0, np.sqrt(3), 0],   
        [-1.0, 0, ca],
    ]) * a
    
    motif = np.array([
        [  0,    0,    0],
        [1/2,  1/2,    0],
        [  0,  5/6,  1/2],
        [1/2,  1/3,  1/2]
    ])

    ncell2 = ncell.copy()
    for i in np.arange(3):
        ncell2[i] = ncell[ np.mod(i+2, 3) ]
   
    atoms = create_supercell_fast(latt, motif, ncell2)
    atoms = vf.make_SFP_xy(atoms, i1=1)
    atoms = vf.make_a3_ortho(atoms)
    

    if bp == 33:
        atoms.positions = atoms.positions + np.array([0, 0, 0.1])*a
        atoms.wrap()

    atoms.pos_a0 = a 
    return atoms
    #vf.my_write_vasp(atoms, filename='POSCAR', vasp5=True)

##### End of code taken from


##### This part is taken from 
# https://github.com/BinglunYin/myalloy_package/blob/master/myvasp/vasp_create_fcc.py

def vasp_create_fcc_100(a, ncell, bp=33):
    """Build an FCC (100) supercell via lattice + motif.

    Uses a cubic lattice (a, a, a) with a 4-atom FCC motif at
    (0,0,0), (0.5,0.5,0), (0.5,0,0.5), (0,0.5,0.5).
    The ``bp`` flag shifts atoms in +z by 0.1*a.

    Args:
        a (float): FCC lattice constant (Å).
        ncell (np.ndarray): (nx, ny, nz) tiling factors.
        bp (int): 33 shifts +z, 0 does not.

    Returns:
        ase.Atoms: FCC (100) supercell.

    Note:
        Depends on the optional ``myvasp`` companion package.
    """
    print('==> create fcc 100 plane: ')
    print(a, ncell, bp)

    latt = np.array([
        [1.0, 0, 0],
        [0, 1.0, 0],   
        [0, 0, 1.0],
    ]) * a
   
    motif = np.array([
        [0, 0, 0],
        [0.5,  0.5,  0],
        [0.5,  0,  0.5],
        [0,  0.5,  0.5],
    ])

    atoms = create_supercell_fast(latt, motif, ncell)
    atoms.pos_a0 = a 

    if bp == 33:
        atoms.positions = atoms.positions + np.array([0, 0, 0.1]) *a 
        atoms.wrap()

    return atoms
    #vf.my_write_vasp(atoms, filename='POSCAR', vasp5=True)

##### End of code taken from