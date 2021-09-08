def Help_GenomeGenerator(GenesDF, GenomeSize, GCcont) -> str :
    '''
    Constructs whole genome with interspersed genes.
    '''
    import numpy as np
    from ..random import pick_sample
    # combining promoter and ORF
    Genes = [''.join([myProm, myORF]) for myProm, myORF in zip(GenesDF['Promoter'].values,GenesDF['ORF'].values)]
    Genes_List = [Convert(Gene) for Gene in Genes]
    # generating background genome sequence
    Genome_Bckgd = make_GenomeBckgd(GenomeSize, GCcont)
    # determining position for gene insertion
    Gene_Positions = np.sort(pick_sample(range(len(Genome_Bckgd)),len(Genes)))
    # breaking the background genome in nested lists at gene positions
    Genome_Tmp = make_NestedList(Genome_Bckgd, Gene_Positions)
    # Now inserting the genes
    Gtmp = [np.concatenate([mbg,bed]) for mbg,bed in zip(Genome_Tmp[:-1],Genes_List)]
    Gtmp = np.concatenate([Gtmp,Genome_Tmp[-1]])
    Genome = ''.join([''.join(elm) for elm in Gtmp])

    return Genome

def make_GeneJoiner(HostName, Model, RctID):
    '''
    Determines promoter activity and combines it with enzyme id and ORF.
    Output
        Gene_Info:    dictionary, gene id, gene expression, promoter, ORF
    '''
    from .expression import make_Promoter, Help_PromoterStrength

    Gene_ORF = make_ORF(Model)
    Gene_Promoter = make_Promoter()
    Gene_Activity = Help_PromoterStrength(HostName, Gene_Promoter, Similarity_Thresh=.8)

    Gene_Dict = {'RctID': RctID, 'Expression': Gene_Activity, 'Promoter': Gene_Promoter, 'ORF': Gene_ORF}

    return Gene_Dict

def make_GenomeBckgd(GenomeSize, GCcont):
    '''
    Function for setup of background genome.
    Input
        GenomeSize:      integer, genome nucleotide number
        GCcont:     float, [0..1], GC-content approximate
    Output
        Genome_Bckgd: string
    '''
    from ..random import pick_choice

    SeqNestList = [pick_choice([Letter for Nest in pick_choice([['G','C'],['A','T']], weights=[GCcont, 1-GCcont]) for Letter in Nest]) for _x in range(GenomeSize)]
    Genome_Bckgd = ''.join(Letter for Nest in SeqNestList for Letter in Nest)

    return Genome_Bckgd

def make_ORF(Model):
    '''
    Function to generate a single ORF. Triplet frequencies are taken from E.coli:
    https://openwetware.org/wiki/Escherichia_coli/Codon_usage
    ORF starts always with 'ATG'.
    Input
        CodonTriplets:     dataframe, automatic load, base triplets, ID (e.g. 'Stop', 'Met'), frequency in percent
        Mutant:            class, subfield model.reactions is used to count the enzyme number to derive minimum sequence length
    Output
        Gene_ORF:          string, open reading frame of enzyme
    '''
    from ..random import pick_choice
    import numpy as np
    import pandas as pd

    # coding sequence construction
    # first we determine the minimum coding gene length of nucleotides to distinguish the enzymes in the model
    Enzyme_Number = len(Model.reactions)
    Gen_Minimum = np.ceil(np.log2(Enzyme_Number))

    # we want to represent codon triplicates, we calculate the next highest divisor of three
    Gene_Length = int(np.ceil(Gen_Minimum/3))
    CodonFile = 'CodonTriplets.csv' # TODO: Maybe not found, could be in data folder.
    CodonTriplets = pd.read_csv(CodonFile, delimiter=';', skipinitialspace=True)
    CodonStop = CodonTriplets[['Stop' in s for s in CodonTriplets['Name']]].reset_index()
    CodonCoding = CodonTriplets.drop(CodonStop.index).reset_index()
    Gene_ORF = [pick_choice(CodonCoding['Triplet'], weights=CodonCoding['Percent']) for CodonId in range(Gene_Length)]
    Gene_Stop = pick_choice(CodonStop['Triplet'], weights=CodonStop['Percent'])
    Gene_ORF.append(Gene_Stop)
    Gene_ORF.insert(0,['ATG'])
    Gene_ORF = ''.join([Letter for Nest in Gene_ORF for Letter in Nest])

    return Gene_ORF


def make_NestedList(List,Breaks):
    '''
    Generates a nested list at the break points
    '''
    import numpy as np
    # adding the last element of the list, otherwise elements after the last break disappear
    Breaks = np.concatenate([Breaks,[len(List)]])
    mystart = 0
    List_Nested = list()
    for myend in Breaks:
        List_Nested.append([List[myidx2] for myidx2 in range(mystart,myend)])
        mystart = myend

    return List_Nested

def Convert(string):
    '''
    https://www.geeksforgeeks.org/python-program-convert-string-list/
    '''
    list1=[]
    list1[:0]=string
    return list1

def list_integer(SeqList):
    '''define input values'''
    alphabet = 'ACGT'
    char_to_int = dict((c,i) for i,c in enumerate(alphabet))
    IntegerList = list()
    for mySeq in SeqList:
        # integer encode input data
        integer_encoded = [char_to_int[char] for char in mySeq.upper()]
        IntegerList.append(integer_encoded)
    return IntegerList

def list_onehot(IntegerList):
    OneHotList = list()
    for integer_encoded in IntegerList:
        onehot_encoded = list()
        for value in integer_encoded:
            letter = [0 for _ in range(4)]
            letter[value] = 1
            onehot_encoded.append(letter)
        OneHotList.append(onehot_encoded)
    return OneHotList

