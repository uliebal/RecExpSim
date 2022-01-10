def Help_Test(Host):
    '''
    Testing wether everything works
    '''
    import numpy as np
    from BioLabSim.ModuleMeasureOrganism.Physiology import Express_Max

    try:
        Host.make_TempGrowthExp(np.array([Host.opt_growth_temp]))
        print('Functional: Temperature growth experiment')
    except:
        print('Error in \'make_TempGrowthExp\'')

    try:
        Clone_ID1 = 'Clone_1'
        Promoter1 = 'GCCCATTGACAAGGCTCTCGCGGCCACCTATAATTGCACG'
        Primer1 =   'CGGGTAACTGTTCCGAGAG'
        Tm = 62 # melting temperature
        # cloning:
        Host.make_Cloning(Clone_ID1, Promoter1, Primer1, Tm)
        print('Functional: Cloning experiment')
    except:
        print('Error in \'make_Cloning\'')

    try:
        Host.measure_PromoterStrength('Clone_1')
        print('Functional: Promoter measurement')
    except:
        print('Error in \'measure_PromoterStrength\'')

    try:
        Host.measure_ProductionExperiment('Clone_1', Host.opt_growth_temp, 1, Host.max_biomass)
        achieveExpRate = Express_Max(Host)

        if Host.var_Library['Clone_1']['Expression_Rate'] == achieveExpRate:
            print('Functional: Expression experiment')
        else:
            print('Warning: Expression experiment with different results')
    except:
        print('Error in \'make_Cloning\'')
###################################################
#
###################################################
def Help_Progressbar(n, loading_time, add):
    '''function for display of a loading bar, n: width of loading bar'''
    import sys
    import time

    loading = '.' * n
    for i in range(n+1):
        # this loop replaces each dot with a hash!
        print('\r%s progress{}: %3d percent'.format(add) % (loading, i*100/n), end='')
        loading = loading[:i] + '#' + loading[i+1:]
        time.sleep(loading_time)
    sys.stdout.write("\n")
###################################################
#
###################################################
def Error_Resources():
    print('Not enough resources available.')
###################################################
#
###################################################
def Help_SwitchComplementary(argument):
    switcher = {
        'T': 'A',
        'A': 'T',
        'C': 'G',
        'G': 'C'
    }
    return switcher.get(argument)
###################################################
#
###################################################
def Sequence_ReferenceDistance(SeqObj, RefSeq=None):
    '''Returns the genetic sequence distance to a reference sequence.
    Input:
           SeqDF: list, the sequence in conventional letter format
    Output:
           SequenceDistance: float, genetic distances as determined from the sum of difference in bases divided by total base number, i.e. max difference is 1, identical sequence =0
    '''
    import numpy as np

    if RefSeq != None:
        RefSeq = SeqObj[0]
    else:
        RefSeq = 'GCCCATTGACAAGGCTCTCGCGGCCAGGTATAATTGCACG'

    Num_Samp = len(SeqObj)
    SequenceDistance = np.sum([int(seq1!=seq2) for seq1,seq2 in zip(RefSeq, SeqObj)], dtype='float')/len(SeqObj)

    return SequenceDistance
###################################################
#
###################################################
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
###################################################
#
###################################################
