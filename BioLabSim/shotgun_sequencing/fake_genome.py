
from typing import List, Optional
from bisect import bisect_left

from numpy.random import Generator as RandomGenerator, PCG64

from .common import Base, Sequence



class SimplifiedGenome :

    # Base sequence on the 5' to 3' strand.
    # TODO: Is this name correct?
    template_strand : Sequence

    def __init__ ( self ) :
        self.template_strand = []



class RandomGenome (SimplifiedGenome) :

    def __init__ ( self, gc_content:float = 0.5, num_bp:int = 100, seed:Optional[int] = None ) :
        """
        Arguments
            gc_content  Probability between 0 and 1 (inclusive) of having G or C on each base.
            num_bp      Number of base pairs.
            seed        Specify a randomization seed for pseudorandom generation.
        """
        SimplifiedGenome.__init__(self)

        # Initialize the seed for random generation.
        randgen = RandomGenerator(PCG64(seed))

        # Calculate the chance of getting each base. Use cummulative probabilities.
        rnd_base:List[Base] = [ 'A', 'T', 'C', 'G' ]
        rnd_ceil:List[float] = [ 0, 0, 0, 0 ]
        rnd_ceil[0] = (1 - gc_content) * 0.5
        rnd_ceil[1] = (1 - gc_content) * 0.5 + rnd_ceil[0]
        rnd_ceil[2] = gc_content * 0.5 + rnd_ceil[1]
        rnd_ceil[3] = gc_content * 0.5 + rnd_ceil[2] # Should always be 1.

        # Generate the strand following argument rules.
        self.template_strand = [
            rnd_base[ bisect_left( rnd_ceil, randgen.random() ) ]
            for _ in range(num_bp)
        ]






