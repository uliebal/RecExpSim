def measure_BaseCompare(Seq1, Seq2):
    '''
    Comparison of two sequences from start to end. returns positions of base differences.
    '''
    SeqDiff = [[count, Pos] for count, Pos in enumerate(zip(Seq1, Seq2)) if Pos[0] != Pos[1]]
    
    return SeqDiff