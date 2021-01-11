"""
This file contains helper methods to visualize data structures and aid in debugging.
"""

from math import floor, ceil
from typing import Dict, List

import numpy as np

from ..common import Scaffold
from .assembly import LocalizedSequence, to_consensus
from ..common import Sequence



HeatPoint = Dict[str,int]



def print_scaffold ( scaffold:Scaffold ) -> None :
    """
    Print the scaffold content as FASTQ.
    """
    print(scaffold.r1_seqrecord.format("fastq"))
    if scaffold.r2_seqrecord is not None :
        print(scaffold.r2_seqrecord.format("fastq"))



def pretty_print_scaffold ( scaffold:Scaffold ) -> None :
    """
    Visually print a scaffold.
    """
    r2_exists = scaffold.r2_seqrecord != None

    # Calculate the gap size.
    gap = 0
    if scaffold.read_method == 'single-read' :
        gap = scaffold.expected_len - len(scaffold.r1_seqrecord.seq)
    elif scaffold.read_method == 'paired-end' :
        gap = scaffold.expected_len - len(scaffold.r1_seqrecord.seq) - len(scaffold.r2_seqrecord)

    print( "| {r1s}~{lspad}(ca.{gap}){rspad}~{r2s}\n| {r1q}~{lspad}~~~~{cpad}{rspad}~~{r2q}".format(
        r1s="".join(scaffold.r1_seqrecord.seq),
        r1q="".join(map(str,[ min(9,floor(q*10)) for q in scaffold.r1_quality ])),
        r2s="".join(reversed(scaffold.r2_seqrecord)) if r2_exists else "",
        r2q="".join(map(str,[ min(9,floor(q*10)) for q in reversed(scaffold.r2_quality) ])) if r2_exists else "",
        lspad="~" * ceil(gap / 2),
        rspad="~" * floor(gap / 2),
        gap=gap,
        cpad="~" * len(str(gap))
    ))



def numtochar ( num:int ) -> str :
    """
    Convert a number to a visual character.
    """
    if num == 0 :
        return "·"
    elif num >= 1 and num <= 9 :
        return str(num)
    elif num > 9 :
        return "#"
    return "?"



def print_assembly_evaluation( loc_sequences:List[LocalizedSequence], genome:Sequence ) -> None :
    """
    Visual view of the evaluation between the placement of multiple localized sequences and the
    real genome. It help in understanding for assembly may possibly work, though high-quality
    assemblers will not be working with lists of sequences, as they try to form a finalized sequence
    with consensus.
    """

    # Build a heatmap of all the sequences overlapped on themself. Use it to align the bundle of
    # sequences with the genome.
    seq_starts = np.array([ locseq.locus for locseq in loc_sequences ])
    seq_ends = np.array([ locseq.locus + len(locseq.sequence) for locseq in loc_sequences ]) # Excluding end.
    bundle_start = np.amin(seq_starts) # locus where the bundle of localized sequences starts
    bundle_end = np.amax(seq_ends)
    bundle_len = bundle_end - bundle_start
    heatmap:List[HeatPoint] = [ {'A':0,'C':0,'G':0,'T':0,'-':0} for _ in range(bundle_len) ]

    # Fill the heatmap with each base that was read.
    for locseq in loc_sequences :
        for s_i in range(len(locseq.sequence)) : # s_i = index on a sequence
            base = locseq.sequence[s_i]
            heatmap[ locseq.locus + s_i - bundle_start ][base] += 1

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
                hm_total = heatmap[b_i]['A'] + heatmap[b_i]['C'] + heatmap[b_i]['G'] + heatmap[b_i]['T'] + heatmap[b_i]['-']
                score += ( heatmap[b_i][gen_base] / hm_total ) if hm_total > 0 else 0

        if score > best_score :
            best_score = score
            best_genome_start = genome_start

    # Calculate the matching score of each singular sequence.
    seq_scores = [ 0 for _ in loc_sequences ]
    for l_i, locseq in enumerate(loc_sequences) :
        score = 0
        for s_i in range(len(locseq.sequence)) :
            g_i = locseq.locus + s_i - best_genome_start
            if g_i in range(genome_len) :
                if locseq.sequence[s_i] == genome[g_i] :
                    score += 1
        seq_scores[l_i] = score

    # Get the consensus sequence.
    consensus = to_consensus(loc_sequences)

    # Display all the sequences in alignment.
    window_start = min( bundle_start, best_genome_start )
    window_end = max( bundle_end, best_genome_start + genome_len )
    print( "-" * (13 + window_end - window_start) )
    for base in ['-','A','C','G','T'] :
        print("[heatmap:{}] {}{}{}".format(
            base,
            numtochar(0) * (bundle_start - window_start),
            "".join([ numtochar(hp[base]) for hp in heatmap ]),
            numtochar(0) * (window_end - bundle_end),
        ))
    print(" [coverage] {}{}{}".format(
        numtochar(0) * (bundle_start - window_start),
        "".join([ numtochar(hp['A']+hp['T']+hp['C']+hp['G']+hp['-']) for hp in heatmap ]),
        numtochar(0) * (window_end - bundle_end),
    ))
    print( "-" * (13 + window_end - window_start) )
    for l_i, locseq in enumerate(loc_sequences) :
        print("  [seq:{:03}] {}{}{} ({:3.0f}%)".format(
            l_i,
            " " * (locseq.locus - window_start),
            "".join([ base for base in locseq.sequence ]),
            " " * (window_end - locseq.locus - len(locseq.sequence)),
            (seq_scores[l_i] / len(locseq.sequence)) * 100
        ))
    print( "-" * (13 + window_end - window_start) )
    print("[consensus] {}{}".format(
        " " * (consensus.locus - window_start),
        "".join([ base for base in str(consensus.sequence) ])
    ))
    print( "-" * (13 + window_end - window_start) )
    print(" [real-gen] {}{}".format(
        " " * (best_genome_start - window_start),
        "".join([ base for base in genome ])
    ))
    print( "-" * (13 + window_end - window_start) )
    print( "Final Score: {:.3} ({:.1f}%)".format(
        best_score,
        (best_score / genome_len) * 100
    ) )
