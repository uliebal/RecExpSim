def make_MutActProm(Genome, GenesDF, NumberEnzymes=3, Target='-10', NumberMutations=3):
    '''
    Add mutations to the promoter of an active enzyme and returns the genome.
    '''
    import random
    FluxActive = GenesDF[GenesDF['Fluxes']!=0].index.values
    MutateEnzyme = random.sample(list(FluxActive),NumberEnzymes)
    
    for i1 in range(NumberEnzymes):
        RefProm = GenesDF['Promoter'].iloc[MutateEnzyme[i1]]
        # mutations in -10box
        if Target == '-10':
            # extracting reference -10 box sequence
            B_i, B_s = -12, -7
            Extract = RefProm[B_i:B_s]
            # finding positions to mutate
            Mutate_Pos = random.sample(range(5), NumberMutations)
            # generating new sequence with the remaining nucleotides at each position
            [Base for Base in mylist if Base is not Extract[Pos]] for Pos in Mutate_Pos