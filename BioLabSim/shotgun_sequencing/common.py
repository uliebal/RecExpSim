
from typing import List, Literal, NamedTuple



# One of the 4 bases.
Base = Literal['A','T','C','G']

# A sequence of bases.
Sequence = List[Base]

# A sequence with a starting position.
class LocalizedSequence (NamedTuple) :
    sequence : Sequence
    locus : int

