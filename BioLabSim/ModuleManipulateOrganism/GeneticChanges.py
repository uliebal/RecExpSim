def Help_MutActProm(Genome, GenesDF, NumberEnzymes=3, Target='-10', NumberMutations=2):
    '''
    Add mutations to the promoter of an active enzyme and returns the genome.
    '''
    import random
    
    FluxActive = GenesDF[GenesDF['Fluxes']!=0].index.values
    MutateEnzyme = random.sample(list(FluxActive),NumberEnzymes)
#     print('Mutated Enzymes: {}'.format(MutateEnzyme))
    Genome_Mutated = Genome
    GenesDF_Mutated = GenesDF.copy()
    GenesDF_Mutated.drop(columns=['Fluxes','Expr2Flux'], inplace=True)
    
    for i1 in range(NumberEnzymes):
        RefProm = GenesDF['Promoter'].iloc[MutateEnzyme[i1]]
#         RefORF = GenesDF['ORF'].iloc[MutateEnzyme[i1]]
        # mutations in -10box
        if Target == '-10':
            # extracting reference -10 box sequence
            B_i, B_s = -13, -6
            RefTar = RefProm[B_i:B_s]
        elif Target == '-35':
            # -35 box region
            B_i, B_s = -37, -30
            RefTar = RefProm[B_i:B_s]        
        else: 
            # whole promoter region
            RefTar = RefProm
        
        MutTar = make_Mutate(RefTar, NumberMutations)
        MutProm = RefProm.replace(RefTar, MutTar)
        Genome_Mutated = Genome_Mutated.replace(RefProm, MutProm)
        GenesDF_Mutated.loc[MutateEnzyme[i1], 'Promoter'] = MutProm

    return Genome_Mutated, GenesDF_Mutated
            
            
def make_Mutate(Sequence, NumberMutations=2):
    '''
    Insert mutations in a given sequence
    '''
    import random

    Bases = ['A','C','G','T']
    
    MutTar = list(Sequence)
    # finding positions to mutate
    Mutate_Pos = random.sample(range(len(Sequence)), NumberMutations)
    # generating new sequence with the remaining nucleotides at each position
    for NuclPos in Mutate_Pos:
        MutTar[NuclPos] = random.sample([Base for Base in Bases if Base is not Sequence[NuclPos]], 1)[0]

    return ''.join(MutTar)
