
from typing import List

import numpy as np

from .common import Sequence, LocalizedSequence



OverlappingSequences = List[LocalizedSequence]



class SimpleContigAssembler :
    """
    This DNA fragment assembler will try to match the extremities of the received fragmented
    sequences. It is not optimized at all and uses some simples metrics and minimum matching
    length and difference rate between the entire matching. Disregards quality.
    """

    # Minimum number of matching pairs to consider a join.
    min_match: int

    # Permitted error rate while matching pairs. Number between 0 and 1.
    error_threshold : float



    def __init__ ( self, min_match:int = 8, error_threshold:float = 0.1 ) :
        self.min_match = min_match
        self.error_threshold = error_threshold



    # DEFUNCT - needs to be completed.
    def apply ( self, frames:List[Sequence] ) -> List[OverlappingSequences] :
        """
        Calculate a [L,R] matrix of matching score between every fragment.
            L----------
                   ++++          (score of 4)
                   R----------
        Matching will continue until its not possible anymore or right before the error threshold
        has been reached.
        Very slow, but not greedy and its consistent.
        """

        # Convert the frame sequences into numpy arrays for slight fast runs and simpler comparisons.
        num_f = len(frames)
        f_lens = np.array([ len(frame) for frame in frames ])
        max_f_len = np.amax(f_lens)
        seq_mat = np.zeros( (num_f,max_f_len), dtype=int )
        # TODO: Complete this.

        score = np.zeros( (numf,numf), dtype=int )
        for a in range(numf) :
            for b in range(numf) :
                if a != b :
                    len_a = len(frames[a])
                    len_b = len(frames[b])

                    # Compare the B string, starting at the leftmost position until the rightmost.
                    for pad_b in range( 0, len_a - self.min_match + 1 ) : # TODO: Maybe -1
                        pass # TODO: Complete this.

        # TODO: Complete this.





