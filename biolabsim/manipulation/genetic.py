
def Help_MutActProm(HostName, Genome, GenesDF, NumberEnzymes=3, Target='-10', NumberMutations=2):
    '''
    Add mutations to the promoter of an active enzyme and returns the genome.
    '''
    from ..random import pick_sample
    from ..simulation.expression import Help_PromoterStrength

    FluxActive = GenesDF[GenesDF['Fluxes']!=0].index.values
    MutateEnzyme = pick_sample(list(FluxActive),NumberEnzymes)
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
        Gene_Activity = Help_PromoterStrength(HostName, MutProm, Similarity_Thresh=.8)
        Genome_Mutated = Genome_Mutated.replace(RefProm, MutProm)
        GenesDF_Mutated.loc[MutateEnzyme[i1], 'Expression'] = Gene_Activity
        GenesDF_Mutated.loc[MutateEnzyme[i1], 'Promoter'] = MutProm

    return Genome_Mutated, GenesDF_Mutated


def make_Mutate(Sequence, NumberMutations=2):
    '''
    Insert mutations in a given sequence
    '''
    from ..random import pick_sample

    Bases = ['A','C','G','T']

    MutTar = list(Sequence)
    # finding positions to mutate
    Mutate_Pos = pick_sample(range(len(Sequence)), NumberMutations)
    # generating new sequence with the remaining nucleotides at each position
    for NuclPos in Mutate_Pos:
        MutTar[NuclPos] = pick_sample([Base for Base in Bases if Base is not Sequence[NuclPos]], 1)[0]

    return ''.join(MutTar)

def Help_Cloning(Host:'Host', Clone_ID, Promoter, Primer, Tm):
    '''Experiment to clone selected promoter. It is displayed whether the experiment was successfull.'''
    import numpy as np
    from ..random import pick_uniform
    from ..auxfun import Help_SwitchComplementary, Sequence_ReferenceDistance


    if Sequence_ReferenceDistance(Promoter) > .4:
        return print('Promoter sequence deviates too much from the given structure.')

    if Host.resources > 0:

        NaConc = 0.1 # 100 mM source: https://www.genelink.com/Literature/ps/R26-6400-MW.pdf (previous 50 mM: https://academic.oup.com/nar/article/18/21/6409/2388653)
        OptLen = Host.opt_primer_len
        AllowDevi = 0.2 # allowed deviation
        Primer_Length = len(Primer)
        Primer_nC = Primer.count('C')
        Primer_nG = Primer.count('G')
        Primer_nA = Primer.count('A')
        Primer_nT = Primer.count('T')
        Primer_GC_content = ((Primer_nC + Primer_nG) / Primer_Length)*100 # unit needs to be percent

        Primer_Tm_1 = 81.5 + 16.6*np.log10(NaConc) + 0.41*Primer_GC_content - 600/Primer_Length # source: https://www.genelink.com/Literature/ps/R26-6400-MW.pdf (previous: https://core.ac.uk/download/pdf/35391868.pdf#page=190)
        Primer_Tm_2 = (Primer_nT + Primer_nA)*2 + (Primer_nG + Primer_nC)*4
        # Product_Tm = 0.41*(Primer_GC_content) + 16.6*np.log10(NaConc) - 675/Product_Length
        # Ta_Opt = 0.3*Primer_Tm + 0.7*Product_Tm - 14.9
        # source Product_Tm und Ta: https://academic.oup.com/nar/article/18/21/6409/2388653
        # Product_Length would be the length of the promoter (40)? too small -> negative number comes out for Product_Tm

        error_1 = pick_uniform(-1,1)*0.1*Primer_Tm_1
        error_2 = pick_uniform(-1,1)*0.1*Primer_Tm_2
        Primer_Tm_err_1 = error_1 + Primer_Tm_1
        Primer_Tm_err_2 = error_2 + Primer_Tm_2

        DeviLen = np.absolute(OptLen - Primer_Length)/OptLen
        DeviTm_1 = np.absolute(Primer_Tm_err_1 - Tm)/Primer_Tm_err_1
        DeviTm_2 = np.absolute(Primer_Tm_err_2 - Tm)/Primer_Tm_err_2
        DeviTm = min(DeviTm_1, DeviTm_2)

        #create the complementary sequence of the primer to check for mistakes:
        PrimerComp = ""
        for base in Primer:
            PrimerComp = PrimerComp + Help_SwitchComplementary(base)

        if DeviLen <= AllowDevi and DeviTm <= AllowDevi/2 and Primer_Length <= 30 and PrimerComp == Promoter[:len(Primer)]:
            print('Cloning was successfull.')
            Host.var_Library[Clone_ID] = {}
            Host.var_Library[Clone_ID]['Promoter_Sequence'] = Promoter
            Host.var_Library[Clone_ID]['Promoter_GC-content'] = (Promoter.count('C') + Promoter.count('G')) / len(Promoter)

        else:
            print('Cloning failed')

        Host.resources -= 1

    else:
        print('Not enough resources available.')

