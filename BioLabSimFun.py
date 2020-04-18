# Biotechnology Laboratory Simulator Functions

class Mutant:
    '''The 'Mutant' class stores all information about the organism and the integrated recombinant protein.''' 
    
#     from BioLabSimFun import SequenceRandomizer_Single
    from random import randint
    # random assignment of the production phase, either during growth phase or stationary phase
    __ProdPhase = 'exponential' if randint(0,1)==0 else 'stationary'
    # resources, e.g. money, for conducting tests
    __Resources = 10
    __BiomassMax = None
    
    def __init__(self, Host):
        from random import randint
        self.var_Host = Host
#         self.var_Promoter = []
        self.var_Resources = self._Mutant__Resources
        # optimal growth temperature, randomly assigned
        self.__OptTemp = randint(25,40) # unit: degree celsius
        # maximum biomass concentration, can be adjusted later, now randomly set
        if self.var_Host == 'Ecol':
            self.__BiomassMax = randint(10,100) # unit: in gCDW/l
        elif self.var_Host == 'Pput':
            self.__BiomassMax = randint(60,150) # unit: in gCDW/l
    
    def show_BiotechSetting(self):
        '''Report of all properties defined in the biotech experiment.'''
        self.var_Resources = self._Mutant__Resources
        MyVars = [i for i in list(vars(self).keys()) if 'var_' in i]
        for i in range(len(MyVars)):
            print('{}: {}'.format(MyVars[i].replace('var_',''), getattr(self, MyVars[i])))        
    
    
    def add_Promoter(self, Promoter):
        self.var_Promoter = Promoter
        self.var_GC_content = (self.var_Promoter.count('C') + self.var_Promoter.count('G')) / len(self.var_Promoter)

    def __add_RandomPromoter(self):
        self.var_Promoter = SequenceRandomizer_Single()
        self.var_GCcontent = (self.var_Promoter.count('C') + self.var_Promoter.count('G')) / len(self.var_Promoter)
        
    def Make_MeasurePromoterStrength(self):
        if not hasattr(self, 'var_Promoter'):
            print('Error, not promoter added. Add first a promoter "add_promoter(<Sequence>)".')
            return
        if self._Mutant__Resources > 0:
            self.var_PromoterStrength = Help_Expression(self)
            self._Mutant__Resources -= 1
        else:
            print('Not enough resources available.')
            self.var_PromoterStrength = None
        
    def Make_ProductionExperiment(self, CultTemp):
        if not hasattr(self, 'var_PromoterStrength'):
            print('Error, no promoter available. Add first a promoter and measure promoter strength "Make_MeasurePromoterStrength".')
            return
        
        if self._Mutant__Resources > 0:
            r = Help_GrowthConstant(self, CultTemp)
            GrowthMax = Growth_Maxrate(self, r)
            self._Mutant__Resources -= 1
            self.var_ExpressionRate = round(GrowthMax * self.var_PromoterStrength,2)
        else:
            print('Not enough resources available.')
            self.var_ExpressionRate = None
            
            
    def Make_TempGrowthExp(self, CultTemps, draw_plot=False):
        '''Experiment to determine optimal growth rate. The experiment runs until the maximum biomass is reached.'''
        import numpy as np
        import pandas as pd
        import matplotlib.pyplot as plt
        import pylab as pl
        from IPython import display
        import time 
#        from IPython import embed
#         %matplotlib inline
        
    
        if self._Mutant__Resources > 0:
            
            capacity = self._Mutant__BiomassMax
            # the time of the half maximum population (inflection point) is calculated according to here:
            # https://opentextbc.ca/calculusv2openstax/chapter/the-logistic-equation/
            d_mult = 2 # we multiply the inflection point with 'd_mult' to increase cultivation time
            P0 = 0.1
            
            # determine time vector with maximum length:
            OptTemp = self._Mutant__OptTemp
            Temp_tmax = CultTemps[np.argmax(np.absolute(CultTemps-OptTemp))]
            r_tmax = Help_GrowthConstant(self, Temp_tmax)
            duration_tmax = d_mult * 1/r_tmax * np.log((capacity - P0)/P0)
            t_max = np.arange(duration_tmax)
            
            # create an empty DataFrame with t_max as first column
            col = []
            col.append("time")
            for i in range (len(CultTemps)):
                col.append(f'biomass {CultTemps[i]}')

            df = pd.DataFrame(np.empty(shape=(len(t_max), len(CultTemps)+1), dtype=float), columns = col)
            df[:len(t_max)] = np.nan 
            new_df = pd.DataFrame({'time': t_max})
            df.update(new_df)
            
            #computing of biomass data and updating of DataFrame
            for i in range (len(CultTemps)):
                r = Help_GrowthConstant(self, CultTemps[i])
                duration = d_mult * 1/r * np.log((capacity - P0)/P0)
                t = np.arange(duration)
                exp_TempGrowthExp = capacity / (1 + (capacity-P0) / P0 * np.exp(-r * t))
                new_df = pd.DataFrame({f'biomass {CultTemps[i]}': exp_TempGrowthExp})
                df.update(new_df)
            
                self._Mutant__Resources -= 1
            
            excel_writer = pd.ExcelWriter('Strain_characerization.xlsx') # Export DataFrame to Excel
            df.to_excel(excel_writer, sheet_name='different temp')
            excel_writer.close()
            
        else:
            print('Not enough resources available.')
            return
        
        if draw_plot:
            WaitingTime = 30/r
            for i in range (len(t)):
                pl.clf()
                pl.figure(figsize = (5,3), dpi = 120)
                pl.xlim(0, t[i])
                pl.ylim(0, 1.1*np.max(exp_TempGrowthExp))
                pl.xlabel("time [h]")
                pl.ylabel("biomass concentration [g/L]")
                pl.plot(t[0:i], exp_TempGrowthExp[0:i], linestyle = '--') 
                display.display(pl.gcf())
                display.clear_output(wait=True)
                #pl.pause(1)
                time.sleep(WaitingTime/n) # total waiting time independent of n  
        
        #Help_ExportToExcel(self, 'Strain_characterization', 'different Temp',
                  #t, exp_TempGrowthExp, 'time', f'biomass {CultTemp}') 
        
    
    def Make_Cloning(self, Primer, Tm):
        '''Experiment to clone selected promoter. The efficiency of the experiment is displayed.'''
        import numpy as np
        import random
        
        if self._Mutant__Resources > 0:
            
            NaConc = 5e-02 # 50 mM laut: https://academic.oup.com/nar/article/18/21/6409/2388653
            Primer_Length = len(Primer)
            Primer_nC = Primer.count('C')
            Primer_nG = Primer.count('G')
            Primer_GC_content = (Primer_nC + Primer_nG) / Primer_Length
            
            Primer_Tm = 81.5 + 16.6*np.log10(NaConc) + 0.41*Primer_GC_content - 675/Primer_Length # Quelle: https://core.ac.uk/download/pdf/35391868.pdf#page=190
            # Product_Tm = 0.41*(Primer_GC_content) + 16.6*np.log10(NaConc) - 675/Product_Length
            # Ta_Opt = 0.3*Primer_Tm + 0.7*Product_Tm - 14.9
            # Quelle Product_Tm und Ta: https://academic.oup.com/nar/article/18/21/6409/2388653
            # Product_Length wäre die Länge des Promoters (40)? zu klein -> negative Zahl kommt raus für Product_Tm
            
            error = random.uniform(-1,1)*0.1*Primer_Tm
            Primer_Tm_err = error + Primer_Tm
            
            exp_Cloning = (1 - np.absolute(Primer_Tm_err - Tm)/Primer_Tm_err) * 100
            print(f'The efficiency of cloning is {exp_Cloning.round(2)} %.')
        
        else:
            print('Not enough resources available.')
            
            
   
        
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

def Help_Expression(Mutant, Predict_File=None, Similarity_Thresh=.4):
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
    
    if Sequence_ReferenceDistance(Mutant.var_Promoter) > Similarity_Thresh:
        Expression = 0
    else:
        if Predict_File!=None:
            Regressor_File = Predict_File
        else:    
            Data_Folder = 'ExpressionPredictor'
            if Mutant.var_Host == 'Ecol':
                Regressor_File = os.path.join(Data_Folder,'Ecol-Promoter-predictor.pkl')
                Add_Params = os.path.join(Data_Folder,'Ecol-Promoter-AddParams.pkl')
                Scaler_DictName = 'Ecol Promoter Activity_Scaler'
            elif Mutant.var_Host == 'Pput':
                Regressor_File = os.path.join(Data_Folder,'Ptai-Promoter-predictor.pkl')
                Add_Params = os.path.join(Data_Folder,'Ptai-Promoter-AddParams.pkl')
                Scaler_DictName = 'Ptai Promoter Activity_Scaler'
            else:
                print('Non-recognized host name. Rename host to either "Ecol" or "Pput."')

        Predictor = joblib.load(Regressor_File)
        Params = pickle.load(open(Add_Params, 'rb'))
        Positions_removed = Params['Positions_removed']
        Expr_Scaler = Params[Scaler_DictName]

        X_Test = np.array(list_onehot(np.delete(list_integer(Mutant.var_Promoter),Positions_removed, axis=0))).reshape(1,-1)  
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

def Help_GrowthConstant(Mutant, CultTemp, var=5):
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


def Help_ExportToExcel(Mutant, FileName, sheet_name, x_values, y_values,
                       x_name, y_name):
    '''Function that exports data from an experiment to an Excel file
    Input:
            FileName: string, name of the file
            sheet_name: string, name of the excel sheet
            x_values: array, x-Werte
            y_values: array, y-Werte
            x_name: string, name of the column with the x-values
            y_name: string, name of the column with the y-values
    Output:
            generated Excel file
    '''
    from pathlib import Path
    import pandas as pd
    
    fname = Path(FileName + '.xlsx')
    if fname.is_file(): # Does the file already exists?
        df = pd.read_excel(FileName + '.xlsx', index_col = 0) # read in and add data
        df[y_name] = y_values
    else:
        df = pd.DataFrame({x_name: x_values,
                            y_name: y_values},
                                columns = [x_name, y_name]) # otherwise create a new DataFrame

    excel_writer = pd.ExcelWriter(FileName + '.xlsx') # Export DataFrame to Excel
    df.to_excel(excel_writer, sheet_name=sheet_name)
    excel_writer.close()