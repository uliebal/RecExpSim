
from typing import List, Optional, NamedTuple

from numpy.random import Generator as RandomGenerator, PCG64

from .fake_genome import SimplifiedGenome
from .common import Sequence



class Frame (NamedTuple) :
    """
    A sequence of bases, with added information about read quality.
    """

    sequence : Sequence
    quality : List[float] # Quality of each base from 0 (worst) to 1 (best).

    def __len__ (self) :
        return len(self.sequence)



class Sequencer :
    """
    TODO: Also to note, shotgun sequencing is more alike taking N sequences and breaking them apart.
    This current method of picking sub-sequences from an unlimited pool is most probably very
    different in reality. It still has many of the mechanisms that can be used in the good version.
    """

    # Average number of base pairs each resulting OpenReadFrame will have.
    # TODO: See how frame length behaves in reality.
    orf_length_mean : int

    # Standard deviation in the length of resulting OpenReadFrames.
    # TODO: See how frame length behaves in reality.
    orf_length_sd : float

    # Number of frames that will be obtained.
    # TODO: See how this behaves in reality.
    num_orfs : int

    # Probability that a base read is wrong.
    # If this occurs, the base that was read will be randomly and equally selected (could even be
    # the correct one).
    # TODO: See how this behaves in reality.
    read_error_prob : float



    def __init__ (
        self,
        orf_length_mean:int = 15, orf_length_sd:float = 3.0, num_orfs:int = 20,
        read_error_prob:float = 0.01
    ) :
        self.orf_length_mean = orf_length_mean
        self.orf_length_sd = orf_length_sd
        self.num_orfs = num_orfs
        self.read_error_prob = read_error_prob

        var_a:str = "ATC"
        var_b:float = 25.74
        var_c = var_a.genome / var_b




    def apply ( self, genome:SimplifiedGenome, seed:Optional[int] = None ) -> List[Frame] :


        # Initialize the seed for random generation.
        randgen = RandomGenerator(PCG64(seed))

        # Go for the approach of selecting randomly a start point and length, then check if it is
        # possible. Repeat if not.
        out_frames:List[Frame] = []
        while len(out_frames) < self.num_orfs :
            start = randgen.integers( low=0, high=len(genome.template_strand) ) # [low,high)
            length = int(randgen.normal( self.orf_length_mean, self.orf_length_sd ))

            # Repeat this sampling if it surpassed the possible range.
            if not ( start + length <= len(genome.template_strand) ) :
                continue

            # Extract the slice and then apply a read error on each base.
            cur_sequence = genome.template_strand[ start : (start+length) ]
            cur_quality = [ 1 for _ in range(len(cur_sequence)) ]
            for i in range(len(cur_sequence)) :
                if randgen.random() < self.read_error_prob :
                    err_base = ['A','T','C','G'][ randgen.integers(0,4) ] # Select a base randomly.
                    cur_sequence[i] = err_base
                    cur_quality[i] = 0.25

            cur_frame = Frame( sequence=cur_sequence, quality=cur_quality )
            out_frames.append(cur_frame)

        return out_frames

