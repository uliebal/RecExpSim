
from typing import List, Literal, NamedTuple, Optional

from Bio.Seq import Seq as BioSeq


# One of the simple 4 nucleic bases, plus the unknown base.
Base = Literal['A','C','G','T']


# Use the same sequence as `biopython` does.
Sequence = BioSeq


# Well distinguished promoter sites.
PromoterSite = Literal[ '-10', '-35' ]


# Read methods usable by the sequencer.
ReadMethod = Literal[ 'single-read', 'paired-end' ]

