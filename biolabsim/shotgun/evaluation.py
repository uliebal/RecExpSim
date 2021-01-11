
from typing import Dict, List

import numpy as np

from .datatype import BaseProbDict, LocalizedSequence, EstimatedSequence
from ..common import Sequence, Base



def evaluate_sequence ( estimated_sequence:EstimatedSequence, genome:Sequence ) -> float :
    """
    Compares an estimated sequence with the real genome and returns the true positive rate.
    Return value is a number between 0 (worst score) and 1 (best score).
    """

    # Build the heatmap by just using the estimation probabilities.
    bundle_start = 0
    bundle_end = len(estimated_sequence)
    bundle_len = bundle_end - bundle_start
    heatmap:List[HeatPoint] = [ {'A':0,'C':0,'G':0,'T':0} for _ in range(bundle_len) ]
    for i in range(bundle_len) :
        for base in ['A','C','G','T'] :
            heatmap[i][base] = estimated_sequence.base_probs[base][i]

    # Now shift the genome around the heatmap and get the best overall match.
    # Score is a simple addition of all ratios of the correct base.
    best_score = 0
    best_genome_start = 0
    genome_len = len(genome)
    for genome_start in range( bundle_start - genome_len, bundle_end ) :
        score = 0
        for g_i in range(genome_len) : # g_i = index on genome sequence
            b_i = genome_start + g_i - bundle_start # b_i = map to an index in the bundle
            if b_i in range(bundle_len) : # only add score if it overlaps with the heatmap, 0 otherwise.
                gen_base = genome[g_i]
                hm_total = heatmap[b_i]['A'] + heatmap[b_i]['C'] + heatmap[b_i]['G'] + heatmap[b_i]['T']
                score += ( heatmap[b_i][gen_base] / hm_total ) if hm_total > 0 else 0

        if score > best_score :
            best_score = score
            best_genome_start = genome_start

    return best_score / genome_len
