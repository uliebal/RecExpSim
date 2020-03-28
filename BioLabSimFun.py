# Biotechnology Laboratory Simulator Functions

class Mutant:
    '''The 'Mutant' class stores all information about the organism and the integrated recombinant protein.''' 
    
#     from BioLabSimFun import SequenceRandomizer_Single
    from random import randint
    # optimal growth temperature, randomly assigned
    __OptTemp = randint(25,40) # unit: degree celsius
    # random assignment of the production phase, either during growth phase or stationary phase
    __ProdPhase = 'exponential' if randint(0,1)==0 else 'stationary'
    # resources, e.g. money, for conducting tests
    __Resources = 3
    __BiomassMax = None
    
    def __init__(self, Host):
        from random import randint
        self.Host = Host
        self.Promoter = []
        self.Resources = self._Mutant__Resources
        # maximum biomass concentration, can be adjusted later, now randomly set
        if Host == 'Ecol':
            self.__BiomassMax = randint(10,100) # unit: in gCDW/l
        elif Host == 'Pput':
            self.__BiomassMax = randint(60,150) # unit: in gCDW/l

    
    def add_promoter(self, Promoter):
        self.Promoter = Promoter
        self.GC_content = (self.Promoter.count('C') + self.Promoter.count('G')) / len(self.Promoter)

    def __add_random_promoter(self):
        self.Promoter = SequenceRandomizer_Single()
        self.GC_content = (self.Promoter.count('C') + self.Promoter.count('G')) / len(self.Promoter)
        
    def measure_promoter_strength(self):
        if self._Mutant__Resources > 0:
            self.Promoter_Strength = Expression(self)
            self._Mutant__Resources -= 1
        else:
            print('Not enough resources available.')
            self.Promoter_Strength = None
        
    def Production_Experiment(self, CultTemp):
        if self._Mutant__Resources > 0:
            r = Gen_GrowthConstant(self, CultTemp)
            GrowthMax = Growth_Maxrate(self, r)
            self._Mutant__Resources -= 1
            self.ExpressionRate = round(GrowthMax * self.Promoter_Strength,2)
        else:
            print('Not enough resources available.')
            self.ExpressionRate = None
        
# def Cultivation(Mutant, Time):
#     '''
#     The cultivation of the host is managed here. It determines the shape of the growth curve based on the optimal growth temperature. The output is maximum of  active biomass which is used for maximum production rate and the biomass integral which is used for final product titer.
#     Arguments:
#         Mutant: class, contains optimal growth temperature, production phase
#         Time:   float, represents the applied cultivation time
#     Output:
#         Biomass_props: dictionary, subfields
#             Biomass_max: float, represents the maximum of active biomass for which the production rate was maximum
#             Biomass_integral: float, represents the area of the biomass, with which the product titer can be calculated
#     '''
#     from random import random
#     Biomass_max = random(1,10)
#     Biomass_integral = random(50,100)
    
#     Biomass_props = {'Biomass_max': Biomass_max, 'Biomass_integral': Biomass_integral}
#     return Biomass_props

def Expression(Mutant, Predict_File=None, Similarity_Thresh=.4):
    '''Expression of the recombinant protein.
        Arguments:
            Mutant: class, contains optimal growth temperature, production phase
            Predict_File: string, address of regression file
        Output: 
            Expression: float, expression rate
    '''
    import os
    import numpy as np
    import joblib
    import pickle
    
    if Sequence_ReferenceDistance(Mutant.Promoter) > Similarity_Thresh:
        Expression = 0
    else:
        if Predict_File!=None:
            Regressor_File = Predict_File
        else:    
            Data_Folder = 'ExpressionPredictor'
            if Mutant.Host == 'Ecol':
                Regressor_File = os.path.join(Data_Folder,'Ecol-Promoter-predictor.pkl')
                Add_Params = os.path.join(Data_Folder,'Ecol-Promoter-AddParams.pkl')
                Scaler_DictName = 'Ecol Promoter Activity_Scaler'
            elif Mutant.Host == 'Pput':
                Regressor_File = os.path.join(Data_Folder,'Ptai-Promoter-predictor.pkl')
                Add_Params = os.path.join(Data_Folder,'Ptai-Promoter-AddParams.pkl')
                Scaler_DictName = 'Ptai Promoter Activity_Scaler'
            else:
                print('Non-recognized host name. Rename host to either "Ecol" or "Pput."')

        Predictor = joblib.load(Regressor_File)
        Params = pickle.load(open(Add_Params, 'rb'))
        Positions_removed = Params['Positions_removed']
        Expr_Scaler = Params[Scaler_DictName]

        X_Test = np.array(list_onehot(np.delete(list_integer(Mutant.Promoter),Positions_removed, axis=0))).reshape(1,-1)  
        Y_Test_norm = Predictor.predict(X_Test)
        Expression = round(float(Expr_Scaler.inverse_transform(Y_Test_norm)),3)

    return Expression

    
def SequenceRandomizer_Parallel(RefSeq, Base_SequencePosition, n=1000):
    '''
    This function generates random sequence combinations. It takes the reference sequence and changes nucleotides at positions that have been experimentally tested. Only as much nucleotides are changed to remain within a given sequence distance.
    '''
    import numpy as np
    import multiprocessing
    from joblib import Parallel, delayed
    from ExpressionExpert_Functions import SequenceRandomizer_Single
    
    num_cores = multiprocessing.cpu_count()
    use_core = min([num_cores, n])
  
    Result = Parallel(n_jobs=use_core)(delayed(SequenceRandomizer_Single)(RefSeq, Base_SequencePosition) for idx in range(n))

    return Result # Sequence_multiple

def SequenceRandomizer_Single():
    '''Generates a randomized sequence for expression evaluation.
        Argument:
            None
        Output:
            Sequence_Single: string, nucleotide sequence
    '''
    import numpy as np
    import random
    import pandas as pd
    
    Alphabet = ['A','C','G','T']
    RefSeq = 'GCCCATTGACAAGGCTCTCGCGGCCAGGTATAATTGCACG'
    Rand_SeqPos = 22
    Base_Pos = pd.DataFrame(np.ones((Rand_SeqPos,4)), columns=Alphabet, dtype=int)
    Base_Pos.index = [(*range(-31,-12,1)),(*range(-10,-7,1))]

    #   Maximum number of nucleotides that can be changed simultaneously based on the sequence distance cut-off
    Nucleotide_Replace_Numb = len(Base_Pos)
    # Generating the positions with nucleotide randomization
    MySynChange = np.array(random.sample(list(Base_Pos.index), Nucleotide_Replace_Numb))
    # the following dataframe has the information of experimentally tested nucleotides as boolean table
    mytmp = Base_Pos.loc[MySynChange]
    # The boolean table of tested nucleotides is converted into an array containing the explicit nucleotide letters
    myArr = np.tile(Alphabet, (Nucleotide_Replace_Numb,1))
    # following, non-tested nucleotides are deleted
    Pos_Del, Nucl_Del = np.where(mytmp.values == 0)
    if not Pos_Del.any(): # != '':
        # deleting non-tested nucleotides
        myArr[tuple([Pos_Del,Nucl_Del])] = 'X'
        
    # Generating a reference sequence to work with
    TstSeq = list(RefSeq)
    # Changing indices from nucleotide oriented to array oriented
    ArSynChange = MySynChange + len(RefSeq)
    # converting the nucleotide array to a list, so we can delete non-tested nucleotides
    Position_list = myArr.tolist()
    Seq_Base = list()
    for Position in Position_list:
        Seq_Base.append(list(set(Position).difference(set('X'))))

        # randomly choosing a possible nucleotide over the total number of exchange positions
    Replace_Bases = [PosIdx[np.random.randint(len(PosIdx))] for PosIdx in Seq_Base]
    #   Replacing the bases in the reference sequence
    for MutIdx, MutNucl in zip(ArSynChange, Replace_Bases):
        TstSeq[MutIdx] = MutNucl    
    Sequence_Single = ''.join(TstSeq)
    
    return Sequence_Single

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

def Gen_GrowthConstant(Mutant, CultTemp, var=5):
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
    
    OptTemp = Mutant._Mutant__OptTemp
    r_pdf = norm(OptTemp, var)
    growth_rate_const = round(r_pdf.pdf(CultTemp),2) / round(r_pdf.pdf(OptTemp),var)
    
    return growth_rate_const
    

def Growth_Maxrate(Mutant, growth_rate_const):
    '''The function calculates the maximum slope during growth.
        Arguments:
            Mutant: class, contains maximum biomass concentration as carrying capacity
            growth_rate_const: float, maximum growth rate constant
        Output:
            growth_rate_max: float, maximum growth rate
    '''
    capacity = Mutant._Mutant__BiomassMax
    # Equation for calculating the maximum slope
    # https://www.tjmahr.com/anatomy-of-a-logistic-growth-curve/
    GrowthMax = capacity*growth_rate_const/4
    
    return GrowthMax
    
    
# def Meas_ProductionRate(Mutant,GrowthMax, Expression):
#     '''The function calculates the maximum production rate. We assume growth depending production, hence maximum production happens during maximum growth rate. The maximum growth rate is calculated based on logistic growth with cultivation temperature. Measuring is a resource demanding process.
#     Arguments:
#         Mutant._Mutant__resources: int, available resources for measurement
#         GrowthMax: float, maximal growth rate
#         Expression: float, expression rate
#     Output:
#         ExpressionRate: float, maximum production rate on biomass
#         '''
    
#     ExpressionRate = GrowthMax * Expression
    
#     return ExpressionRate

##################################################################
##################################################################
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