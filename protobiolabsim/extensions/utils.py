"""
TODO: Almost all of the methods inside this `utils.py` file can me integrated into module functions.
"""


from typing import Tuple

from Bio.Seq import Seq

from ..random import pick_uniform
from ..common import OperationOutcome
from .records.gene.gene import Gene





def Sequence_ReferenceDistance(SeqObj, RefSeq):
    '''Returns the genetic sequence distance to a reference sequence.
    Input:
           SeqDF: list, the sequence in conventional letter format
    Output:
           SequenceDistance: float, genetic distances as determined from the sum of difference in bases divided by total base number, i.e. max difference is 1, identical sequence =0
    '''
    import numpy as np

    Num_Samp = len(SeqObj)
    SequenceDistance = np.sum([int(seq1!=seq2) for seq1,seq2 in zip(RefSeq, SeqObj)], dtype='float')/len(SeqObj)

    return SequenceDistance





def check_primer_integrity (Promoter:Seq, Primer:Seq, Tm:int, RefPromoter:Seq, OptPrimerLen:int) -> OperationOutcome :
    '''
    Experiment to clone selected promoter. Return whether the experiment was successful.
    '''
    import numpy as np

    if Sequence_ReferenceDistance(Promoter, RefPromoter) > .4:
        return OperationOutcome(False,'Promoter sequence deviates too much from the given structure.')

    NaConc = 0.1 # 100 mM source: https://www.genelink.com/Literature/ps/R26-6400-MW.pdf (previous 50 mM: https://academic.oup.com/nar/article/18/21/6409/2388653)
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

    DeviLen = np.absolute(OptPrimerLen - Primer_Length)/OptPrimerLen
    DeviTm_1 = np.absolute(Primer_Tm_err_1 - Tm)/Primer_Tm_err_1
    DeviTm_2 = np.absolute(Primer_Tm_err_2 - Tm)/Primer_Tm_err_2
    DeviTm = min(DeviTm_1, DeviTm_2)

    #create the complementary sequence of the primer to check for mistakes:
    PrimerComp = Primer.complement()

    if DeviLen <= AllowDevi and DeviTm <= AllowDevi/2 and Primer_Length <= 30 and PrimerComp == Promoter[:len(Primer)]:
        return OperationOutcome(True, 'All good')

    if not DeviLen <= AllowDevi :
        return OperationOutcome(False, 'Primer length deviation too big.')
    if not DeviTm <= AllowDevi/2 :
        return OperationOutcome(False, 'Temperature deviation too big.')
    if not Primer_Length <= 30 :
        return OperationOutcome(False, 'Primer length too big.')
    if not PrimerComp == Promoter[:len(Primer)] :
        return OperationOutcome(False, 'Primer not compatible with promoter.')

    return OperationOutcome(False, 'Cloning failed')



def Help_PromoterStrength(
    PromSequence, RefPromoter:str, Scaler=1, Similarity_Thresh=.4, Regressor_File=None, AddParams_File=None
) -> float :
    '''
    TODO: Documentation out of sync.
    Expression of the recombinant protein.
        Arguments:
            Host:       class, contains optimal growth temperature, production phase
            Sequence:     string, Sequence for which to determine promoter strength
            Scaler:       int, multiplied to the regression result for higher values
            Predict_File: string, address of regression file
        Output:
            Expression: float, expression rate
    '''
    import numpy as np
    import joblib
    import pickle

    if Sequence_ReferenceDistance(PromSequence, RefPromoter) > Similarity_Thresh:
        Expression = 0
    else:
        Predictor = joblib.load(Regressor_File)
        Params = pickle.load(open(AddParams_File, 'rb'))
        Positions_removed = Params['Positions_removed']
#         Expr_Scaler = Params[Scaler_DictName]

        X = np.array(list_onehot(np.delete(list_integer(PromSequence),Positions_removed, axis=0))).reshape(1,-1)
        GC_cont = (PromSequence.count('G') + PromSequence.count('C'))/len(PromSequence)
        X = np.array([np.append(X,GC_cont)])
        Y = Predictor.predict(X)
        Expression = round(float(Y)*Scaler,3)

    return Expression





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








def Help_GrowthConstant(OptTemp, CultTemp, var=5):
    '''Function that generates the growth rate constant. The growth rate constant depends on the optimal growth temperature and the cultivation temperature. It is sampled from a Gaussian distribution with the mean at the optimal temperature and variance 1.
    Arguments:
        Opt_Temp: float, optimum growth temperature, mean of the Gaussian distribution
        Cult_Temp: float, cultivation temperature for which the growth constant is evaluated
        var: float, variance for the width of the Gaussian covering the optimal growth temperature
    Output:
        growth_rate_const: float, constant for use in logistic growth equation
    '''

    import numpy as np
    from scipy.stats import norm

    r_pdf = norm(OptTemp, var)
    # calculation of the growth rate constant, by picking the activity from a normal distribution

    growth_rate_const = r_pdf.pdf(CultTemp) / r_pdf.pdf(OptTemp)

    return growth_rate_const





def Growth_Maxrate(growth_rate_const, Biomass):
    '''
    TODO: Documentation is out of sync.
    The function calculates the maximum slope during growth.
        Arguments:
            Host: class, contains maximum biomass concentration as carrying capacity
            growth_rate_const: float, maximum growth rate constant
        Output:
            growth_rate_max: float, maximum growth rate
    '''
    # biomass checks
#     if Biomass > Mutant._Mutantmax_biomass or not Biomass:
#         print('Error, no biomass was set or unexpected value or the maximum possible biomass was exceeded. Enter a value for the biomass again.')

    # Equation for calculating the maximum slope
    # https://www.tjmahr.com/anatomy-of-a-logistic-growth-curve/
    GrowthMax = Biomass * growth_rate_const / 4

    return GrowthMax

