
from typing import List, NamedTuple, Dict, Any

import numpy as np
from Bio import pairwise2
from Bio.Blast.Record import Alignment

from .datatype import Scaffold, LocalizedSequence, EstimatedSequence, get_consensus_from_overlap, \
    estimate_from_overlap
from ..common import Base, Sequence



#Alignment = Any # TODO: Use biopython Alignment typing.



class PairwiseScore (NamedTuple) :
    i : int
    k : int
    score : float



def match ( sa:Sequence, sb:Sequence ) -> Alignment :
    """ Internal method for pairwise matching. """
    return pairwise2.align.localms( sa, sb, 1, 0, -100, -100, gap_char='-' )



class GreedyContigAssembler :
    """
    This assembler will try to assemble the DNA scaffolds by using greedy pairwise matching.
    It will start with the pair of sequences with highest pairwise matching and from there one
    iteratively add one of the remaining sequences by matching it with the consensus sequence of
    all already clustered sequences.

    Attention: Only reads the R1 sequence. Better assembly is achieved with external programs.
    """

    def __init__ ( self ) :
        pass



    def apply_internal ( self, scaffolds:List[Scaffold] ) -> List[LocalizedSequence] :

        # Calculate best pairwise matching score to decide on the first sequences to match.
        fbest = PairwiseScore( 0, 0, -1 )
        for i in range(len(scaffolds)) :
            for k in range( i+1, len(scaffolds) ) :
                cur_align = match( scaffolds[i].r1_seqrecord.seq, scaffolds[k].r1_seqrecord.seq )[0]
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
        cluster.append( LocalizedSequence( sequence=scaffolds[fbest.i].r1_seqrecord.seq, locus=0 ) )
        #print("First: " + "".join(scaffolds[fbest.i].r1_seqrecord.seq))

        # While sequences are still left unclustered.
        while sum(available) > 0 :
            min_start = min([ locseq.locus for locseq in cluster ])
            consensus = get_consensus_from_overlap(cluster).sequence
            #print("Consensus: " + "".join(consensus))
            best = PairwiseScore( 0, 0, -1 )
            best_align:Any = None # TODO: Use biopython Alignment typing.

            # Get the best matching for the cluster consensus.
            for k in range(len(scaffolds)) :
                if available[k] :
                    cur_align = match( consensus, scaffolds[k].r1_seqrecord.seq )[0]
                    #print("      Try: {} ({}) [{}]".format( "".join(scaffolds[k].r1_seqrecord.seq), cur_align.score, k ))
                    #print("         : {}".format("".join(cur_align.seqA)))
                    #print("         : {}".format("".join(cur_align.seqB)))

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
            if best_align.seqB[0] == '-' : # new sequence starts after cluster.
                best_locus = min_start + best_align.start
            elif best_align.seqA[0] == '-' : # new sequence starts before the cluster
                best_locus = min_start - best_align.start
            else : # new sequence start at the same position as the cluster
                best_locus = min_start

            cluster.append( LocalizedSequence( sequence=scaffolds[best.k].r1_seqrecord.seq, locus=best_locus ))
            available[best.k] = False
            #print("    Added: {} [{}] <{},{}>".format(
            #    "".join(scaffolds[best.k].r1_seqrecord.seq), best.k, min_start, best_locus
            #))

        return cluster



    def apply ( self, scaffolds:List[Scaffold] ) -> EstimatedSequence :
        return estimate_from_overlap( self.apply_internal(scaffolds) )
