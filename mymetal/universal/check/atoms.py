"""
atoms submodule

This submodule provides functions for analyzing atomic structures using the Common Neighbor Analysis (CNA) method. It includes
functions for counting different crystal structures (FCC, HCP, BCC, ICO, OTHER) from OVITO files. These functions are designed to
streamline the analysis of atomic configurations in materials science simulations.

Functions:
    - get_cna_count: Get the count of different crystal structures (FCC, HCP, BCC, ICO, OTHER) from an OVITO file.

"""

from ovito.io import import_file
from ovito.modifiers import CommonNeighborAnalysisModifier
from pathlib import Path

def get_cna_count(path_atoms: Path = Path('./')) -> tuple[int, int, int, int, int]:
        """Get the count of different crystal structures (FCC, HCP, BCC, ICO, OTHER) from an OVITO file.
        
        Args:
            path_atoms (Path): The path to the OVITO file containing atomic data.
        
        Returns:
            tuple[int, int, int, int, int]: A tuple containing the counts of FCC, HCP, BCC, ICO, and OTHER structures, respectively.
        """

        if not path_atoms.is_file():
            raise FileNotFoundError(f"The file '{path_atoms}' does not exist.")

        # Read in OVITO for CNA analysis
        pipeline = import_file(Path(path_atoms))
        cna_modifier = CommonNeighborAnalysisModifier(mode=CommonNeighborAnalysisModifier.Mode.AdaptiveCutoff)
        pipeline.modifiers.append(cna_modifier)
        results = pipeline.compute()

        fcc_count = results.attributes['CommonNeighborAnalysis.counts.FCC']
        hcp_count = results.attributes['CommonNeighborAnalysis.counts.HCP']
        bcc_count = results.attributes['CommonNeighborAnalysis.counts.BCC']
        ico_count = results.attributes['CommonNeighborAnalysis.counts.ICO']
        other_count = results.attributes['CommonNeighborAnalysis.counts.OTHER']

        return fcc_count, hcp_count, bcc_count, ico_count, other_count