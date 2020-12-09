
from typing import TypedDict, List

import numpy as np

from .common import LocalizedSequence
from .fake_genome import SimplifiedGenome



class HeatPoint (TypedDict) :
    A : int = 0
    T : int = 0
    C : int = 0
    G : int = 0



def numtochar ( num:int ) -> str :
    if num == 0 :
        return "·"
    elif num >= 1 and num <= 9 :
        return str(num)
    elif num > 9 :
        return "#"
    return "?"



def print_assembly_evaluation( loc_sequences:List[LocalizedSequence], genome:SimplifiedGenome ) -> None :
    """
    Watch that sequences can be translated in their positions.
    """

    # Build a heatmap of all the sequences overlapped on themself. Use it to align the bundle of
    # sequences with the genome.
    seq_starts = np.array([ locseq.locus for locseq in loc_sequences ])
    seq_ends = np.array([ locseq.locus + len(locseq.sequence) for locseq in loc_sequences ]) # Excluding end.
    bundle_start = np.amin(seq_starts) # locus where the bundle of localized sequences starts
    bundle_end = np.amax(seq_ends)
    bundle_len = bundle_end - bundle_start
    heatmap:List[HeatPoint] = [ HeatPoint(A=0,T=0,G=0,C=0) for _ in range(bundle_len) ]

    # Fill the heatmap with each base that was read.
    for locseq in loc_sequences :
        for s_i in range(len(locseq.sequence)) : # s_i = index on a sequence
            base = locseq.sequence[s_i]
            heatmap[ locseq.locus + s_i - bundle_start ][base] += 1

    # Now shift the genome around the heatmap and get the best overall match.
    # Score is a simple addition of all ratios of the correct base.
    best_score = 0
    best_genome_start = 0
    genome_len = len(genome.template_strand)
    for genome_start in range( bundle_start - genome_len, bundle_end ) :
        score = 0
        for g_i in range(genome_len) : # g_i = index on genome sequence
            b_i = genome_start + g_i - bundle_start # b_i = map to an index in the bundle
            if b_i in range(bundle_len) : # only add score if it overlaps with the heatmap, 0 otherwise.
                gen_base = genome.template_strand[g_i]
                score += (
                    heatmap[b_i][gen_base]
                    / ( heatmap[b_i]['A'] + heatmap[b_i]['T'] + heatmap[b_i]['C'] + heatmap[b_i]['G'] )
                )

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
                if locseq.sequence[s_i] == genome.template_strand[g_i] :
                    score += 1
        seq_scores[l_i] = score

    # Display all the sequences in alignment.
    window_start = min( bundle_start, best_genome_start )
    window_end = max( bundle_end, best_genome_start + genome_len )
    print( "-" * (9 + window_end - window_start) )
    for base in ['A','T','C','G'] :
        print("[hmap:{}] {}{}{}".format(
            base,
            numtochar(0) * (bundle_start - window_start),
            "".join([ numtochar(hp[base]) for hp in heatmap ]),
            numtochar(0) * (window_end - bundle_end),
        ))
    print("[coverg] {}{}{}".format(
        numtochar(0) * (bundle_start - window_start),
        "".join([ numtochar(hp['A']+hp['T']+hp['C']+hp['G']) for hp in heatmap ]),
        numtochar(0) * (window_end - bundle_end),
    ))
    print( "-" * (9 + window_end - window_start) )
    for l_i, locseq in enumerate(loc_sequences) :
        print("[seq:{:02}] {}{}{} ({:3.0f}%)".format(
            l_i,
            " " * (locseq.locus - window_start),
            "".join([ base for base in locseq.sequence ]),
            " " * (window_end - locseq.locus - len(locseq.sequence)),
            (seq_scores[l_i] / len(locseq.sequence)) * 100
        ))
    print( "-" * (9 + window_end - window_start) )
    print("[genome] {}{}".format(
        " " * (best_genome_start - window_start),
        "".join([ base for base in genome.template_strand ])
    ))
    print( "-" * (9 + window_end - window_start) )
    print( "Final Score: {:.3} ({:.1f}%)".format(
        best_score,
        (best_score / genome_len) * 100
    ) )


    # for i in range(len(frames)) :
    #     print("[{:03}] {}   ({})".format( i, "".join(frames[i].sequence), len(frames[i]) ))

