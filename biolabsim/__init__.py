"""
The module biolabsim offers methods to create virtual microorganisms and run simulations on their
metabolic networds, promoter strength and second-generation sequencing.
"""

# The `__init__.py` file defines all public-available methods that we allow an user to import.
# Each subfolder has its own `__init__.py` where all methods are declared that should be available
# outside of that sub-folder.
# When a directory has an `__init__.py` file we are telling Python that it is a library.

from .common import Base, Sequence, PromoterSite, ReadMethod

from .random import set_seed

from .simulation import Strain, HostHasNoStrain, Host, Ecol, FabricatedHost, Pput, Entry, Library, Help_GenomeGenerator, Help_PromoterStrength

from .shotgun import (
    Assembler, GreedyContigAssembler, RandomAssembler, Scaffold, LocalizedSequence,
    EstimatedSequence, get_consensus_from_overlap, estimate_from_overlap, calc_total_score,
    calc_sequence_score, evaluate_sequence, Sequencer, write_scaffolds_to_file,
    print_scaffold_as_fastq, print_scaffold, print_assembly_evaluation, print_estimation_evaluation
)

from .measurement import (measure_EnzymeLevel1) # TODO: Fill in with public methods

# from .manipulation import ( ) # TODO: Fill in with public methods
