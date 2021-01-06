
from typing import List, NamedTuple, Dict, Any

from Bio import pairwise2
import numpy as np

from ..common import Base, Sequence, Scaffold



Alignment = Any # TODO: Use biopython Alignment typing.



# A sequence with a starting position.
class LocalizedSequence (NamedTuple) :
    sequence : Sequence
    locus : int



class PairwiseScore (NamedTuple) :
    i : int
    k : int
    score : float



def match ( sa:Sequence, sb:Sequence ) -> Alignment :
    """ Internal method for pairwise matching. """
    return pairwise2.align.localms( sa, sb, 1, 0, -100, -100, gap_char=['X'] )



def to_consensus ( locseqs:List[LocalizedSequence] ) -> Sequence :
    """
    Output the consensus sequence of multiple localized sequences.
    Returns: ( min starting locus, max ending locus, consensus sequence )
    """

    # Pay attention that indexes start at zero while locus could be negative.
    min_start = min([ locseq.locus for locseq in locseqs ])
    max_end = max([ locseq.locus + len(locseq.sequence) for locseq in locseqs ])
    conseq:List[str] = [ "-" for _ in range(min_start,max_end) ]

    # For each position in the consensus sequence, get the value of each localized sequence in that
    # same position. Could be out of bounds.
    for pos in range( min_start, max_end ) :
        quorum:Dict[Base,int] = {}
        for locseq in locseqs :
            if pos >= locseq.locus and pos < locseq.locus + len(locseq.sequence) : # inside bounds

                # Add this base to the possible bases in the consensus. Initialize if required.
                cur_base = locseq.sequence[ pos - locseq.locus ]
                if cur_base not in quorum :
                    quorum[cur_base] = 0
                quorum[cur_base] += 1

        # Based on the quorum, select the best base for this position.
        max_base = max( quorum, key=quorum.get, default='X' ) # Get Key with maximum value. Get first if equal.
        conseq[ pos - min_start ] = max_base

    return conseq





class GreedyContigAssembler :
    """
    This assembler will try to assemble the DNA scaffolds by using greedy pairwise matching.
    It will start with the pair of sequences with highest pairwise matching and from there one
    iteratively add one of the remaining sequences by matching it with the consensus sequence of
    all already clustered sequences.
    TODO: Only reads the R1 sequence for now.
    """

    def __init__ ( self ) :
        pass



    def apply_internal ( self, scaffolds:List[Scaffold] ) -> List[LocalizedSequence] :

        # Calculate best pairwise matching score to decide on the first sequences to match.
        fbest = PairwiseScore( 0, 0, -1 )
        for i in range(len(scaffolds)) :
            for k in range( i+1, len(scaffolds) ) :
                cur_align = match( scaffolds[i].r1_sequence, scaffolds[k].r1_sequence )[0]
                if cur_align.score > fbest.score :
                    fbest = PairwiseScore( i, 0, cur_align.score ) # 0 because k is not used.

        # To simplify coding (at the expense of an additional loop cycle) we will select the fbest.i
        # as the starting sequence and then recalculate the best pairwise match. It should lead to
        # fbest.k again, but we avoid having to write an extra case for this first match.

        # Keep track of all sequences that are left to cluster.
        available:List[bool] = [ True for _ in range(len(scaffolds)) ]
        cluster:List[LocalizedSequence] = []

        # Add the first sequence to start the cluster.
        available[fbest.i] = False
        cluster.append( LocalizedSequence( sequence=scaffolds[fbest.i].r1_sequence, locus=0 ) )
        print("First: " + "".join(scaffolds[fbest.i].r1_sequence))

        # While sequences are still left unclustered.
        while sum(available) > 0 :
            min_start = min([ locseq.locus for locseq in cluster ])
            consensus = to_consensus(cluster)
            print("Consensus: " + "".join(consensus))
            best = PairwiseScore( 0, 0, -1 )
            best_align:Any = None # TODO: Use biopython Alignment typing.

            # Get the best matching for the cluster consensus.
            for k in range(len(scaffolds)) :
                if available[k] :
                    cur_align = match( consensus, scaffolds[k].r1_sequence )[0]
                    print("      Try: {} ({}) [{}]".format( "".join(scaffolds[k].r1_sequence), cur_align.score, k ))
                    print("         : {}".format("".join(cur_align.seqA)))
                    print("         : {}".format("".join(cur_align.seqB)))

                    if cur_align.score > best.score :
                        best = PairwiseScore( 0, k, cur_align.score ) # 0 because i is not used.
                        best_align = cur_align

            # Somehow break the entire method if this is None, should never happen.
            if best_align is None :
                raise "GreedyContigAssembler is not working as intended."

            # From the alignment it is not known if the sequence A starts before or after the
            # sequence B. Unfortunately the `best_align.start` is always positive. Therefore peek
            # which of the sequences starts with the gap character. The one that does will start
            # after the other sequence. All locus is fixed to the locus of the cluster.
            best_locus = min_start
            if best_align.seqB[0] == 'X' : # new sequence starts after cluster.
                best_locus = min_start + best_align.start
            elif best_align.seqA[0] == 'X' : # new sequence starts before the cluster
                best_locus = min_start - best_align.start
            else : # new sequence start at the same position as the cluster
                best_locus = min_start

            cluster.append( LocalizedSequence( sequence=scaffolds[best.k].r1_sequence, locus=best_locus ))
            available[best.k] = False
            print("    Added: {} [{}] <{},{}>".format(
                "".join(scaffolds[best.k].r1_sequence), best.k, min_start, best_locus
            ))


        return cluster








    # # DEFUNCT - needs to be completed.
    # def apply ( self, frames:List[Sequence] ) -> List[OverlappingSequences] :
    #     """
    #     Calculate a [L,R] matrix of matching score between every fragment.
    #         L----------
    #                ++++          (score of 4)
    #                R----------
    #     Matching will continue until its not possible anymore or right before the error threshold
    #     has been reached.
    #     Very slow, but not greedy and its consistent.
    #     """

    #     # Convert the frame sequences into numpy arrays for slight fast runs and simpler comparisons.
    #     num_f = len(frames)
    #     f_lens = np.array([ len(frame) for frame in frames ])
    #     max_f_len = np.amax(f_lens)
    #     seq_mat = np.zeros( (num_f,max_f_len), dtype=int )
    #     # TODO: Complete this.

    #     score = np.zeros( (numf,numf), dtype=int )
    #     for a in range(numf) :
    #         for b in range(numf) :
    #             if a != b :
    #                 len_a = len(frames[a])
    #                 len_b = len(frames[b])

    #                 # Compare the B string, starting at the leftmost position until the rightmost.
    #                 for pad_b in range( 0, len_a - self.min_match + 1 ) : # TODO: Maybe -1
    #                     pass # TODO: Complete this.

    #     # TODO: Complete this.





