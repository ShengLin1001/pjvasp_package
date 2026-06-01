"""
atoms submodule

This submodule provides functions for analyzing atomic structures using the Common Neighbor Analysis (CNA) method. It includes
functions for counting different crystal structures (FCC, HCP, BCC, ICO, OTHER) from OVITO files. These functions are designed to
streamline the analysis of atomic configurations in materials science simulations.

Functions:
    - get_cna_count: Get the count of different crystal structures (FCC, HCP, BCC, ICO, OTHER) from an OVITO file.

"""

from ovito.modifiers import CommonNeighborAnalysisModifier
from pathlib import Path
from ovito.io.ase import ase_to_ovito
from ase import Atoms
from ovito.pipeline import StaticSource, Pipeline

def get_cna_count(atoms: Atoms = None ) -> tuple[int, int, int, int, int]:
        """Get the count of different crystal structures (FCC, HCP, BCC, ICO, OTHER) from an OVITO file.
        
        Args:
            atoms (Atoms): The ASE Atoms object containing atomic data.
        
        Returns:
            tuple[int, int, int, int, int]: A tuple containing the counts of FCC, HCP, BCC, ICO, and OTHER structures, respectively.
        """

        if atoms is None:
            raise ValueError("The 'atoms' argument cannot be None.")

        # Convert ASE Atoms to OVITO format
        atoms_ovito = ase_to_ovito(atoms)

        # Read in OVITO for CNA analysis
        pipeline = Pipeline(source=StaticSource(data=atoms_ovito))
        cna_modifier = CommonNeighborAnalysisModifier(mode=CommonNeighborAnalysisModifier.Mode.AdaptiveCutoff)
        pipeline.modifiers.append(cna_modifier)
        results = pipeline.compute()

        fcc_count = results.attributes['CommonNeighborAnalysis.counts.FCC']
        hcp_count = results.attributes['CommonNeighborAnalysis.counts.HCP']
        bcc_count = results.attributes['CommonNeighborAnalysis.counts.BCC']
        ico_count = results.attributes['CommonNeighborAnalysis.counts.ICO']
        other_count = results.attributes['CommonNeighborAnalysis.counts.OTHER']

        return fcc_count, hcp_count, bcc_count, ico_count, other_count