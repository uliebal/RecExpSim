
from typing import List, Literal, NamedTuple, Optional


# One of the simple 4 nucleic bases, plus the unknown base. (TODO: X here but confirm)
Base = Literal['A','C','G','T','X']


# A sequence of simple bases.
Sequence = List[Base]


# Well distinguished promoter sites.
PromoterSite = Literal[ '-10', '-35' ]


# Read methods usable by the sequencer.
ReadMethod = Literal[ 'single-read', 'paired-end' ]


class Scaffold (NamedTuple) :
    """
    A sequence of bases, with added information about read quality and relative distances.
    Read quality is expressed in a ratio between 0 (worst) and 1 (best).
    Can either have a single contig or 2 contigs with paired ends.
    The R2 contig has a reversed sequence.
    The expected expected_len includes both contigs and the gap.

        |  R1 contig  |            gap            |  R2 contig  |
        [0123=========] - - - - - - - - - - - - - [=========3210]
         ---> order                                   order <---

    """
    read_method : ReadMethod
    expected_len : float
    r1_sequence : Sequence
    r1_quality : List[float]
    r2_sequence : Optional[Sequence]
    r2_quality : Optional[List[float]]

