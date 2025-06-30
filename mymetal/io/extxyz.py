"""extxyz module

This module provides a function to convert Extended XYZ format files into a list of ASE Atoms objects.
It uses the ASE library's robust reading capabilities to handle the Extended XYZ format, allowing for easy manipulation of atomic structures in simulations.

Functions:
    - extxyz_to_atomlist: Converts an Extended XYZ file to a list of ASE Atoms objects.
"""

from ase.io import read


def extxyz_to_atomlist(file: str = None) -> list:
    """Converts Extended XYZ file to list of ASE Atoms objects.
    
    Args:
        file: Path to Extended XYZ format file. If None, may raise error.
    
    Returns:
        List of ASE Atoms objects, one per frame in the trajectory.

    Notes:
        Uses ase.io.read which is more robust than direct extxyz readers.
        For large files, consider processing frames directly from generator.

    Examples:
        >>> atoms_list = extxyz_to_atomlist('trajectory.xyz')  # Get all frames
        >>> len(atoms_list)  # Number of frames
        10
        
        Alternative low-level implementation:
        ```python
        from io import StringIO
        from ase.io.extxyz import read_extxyz

        with open('movie_CONTCAR.xyz', 'r') as f:
            content = f.read()

        atomlist = read_extxyz(StringIO(content), index = -1)  # cant read all frames directly like ':'
        for i, atoms in enumerate(atomlist):
            print(atoms)
        ```
    """
    atomgenerator = read(file, index=':', format='extxyz')
    atomlist = []
    for atoms in atomgenerator:
        atoms.wrap()
        atomlist.append(atoms)
    return atomlist