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
    atoms = move_atoms(atoms, [0,0,0.1*a], if_scale_position=False)
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
        atoms = move_atoms(atoms, [0,0,0.1*a], if_scale_position=False)
    elif mode == 'narrow':
        atoms = move_atoms(atoms, [0,0,-0.05*a], if_scale_position=False)
    atoms.wrap()
    return None

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
    atoms = move_atoms(atoms, [0,0,0.1*a], if_scale_position=False)
    atoms.wrap()
    return atoms

# todo: BCC (110) Plane

