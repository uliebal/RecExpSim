# Biotechnology Laboratory Simulator Functions
# print(vars(myhost))
 
class Mutant:
    '''The 'Mutant' class stores all information about the organism and the integrated recombinant protein.''' 
    
#     from BioLabSimFun import SequenceRandomizer_Single
    from random import randint
    # random assignment of the production phase, either during growth phase or stationary phase
    __ProdPhase = 'exponential' if randint(0,1)==0 else 'stationary'
    # resources, e.g. money, for conducting tests
    __Resources = 10000
    __BiomassMax = None
    __ExpSucRate = None
    
    def __init__(self, Host):
        from random import randint
        self.var_Host = Host
        self.var_Resources = self._Mutant__Resources
        self.var_Substrate = None
        # Library variable containing details to the different tested mutants
        self.var_Library = {}
        # factor which influences the range of the promoter strength, randomly assigned
        self.__InflProStreng = randint(30,50) # explanation see Plot_ExpressionRate 
        # optimal growth temperature, randomly assigned
        self.__OptTemp = randint(25,40) # unit: degree celsius, source: https://application.wiley-vch.de/books/sample/3527335153_c01.pdf
        # optimal Primer length, randomly assigned
        self.__OptPrLen = randint(16,28) # unit: nt, source: https://link.springer.com/article/10.1007/s10529-013-1249-8
        # maximum biomass concentration, the limits for Ecol were set as shown below and the values for Pput were adjusted according to the ratio of the maximum promoter strengths (0.057/0.04) of the optimal sequences (see expression measurement issue).
        if self.var_Host == 'Ecol':
            self.__BiomassMax = randint(30,100) # unit: in gDCW/l, source (german): https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=2&cad=rja&uact=8&ved=2ahUKEwjzt_aJ9pzpAhWGiqQKHb1jC6MQFjABegQIAhAB&url=https%3A%2F%2Fwww.repo.uni-hannover.de%2Fbitstream%2Fhandle%2F123456789%2F3512%2FDissertation.pdf%3Fsequence%3D1&usg=AOvVaw2XfGH11P9gK2F2B63mY4IM
        elif self.var_Host == 'Pput':
            self.__BiomassMax = randint(45,145) # unit: in gDCW/l, source 1: https://onlinelibrary.wiley.com/doi/pdf/10.1002/bit.25474, source 2: https://link.springer.com/article/10.1385/ABAB:119:1:51

    def BuyEquipment(self, EquipInvest=None):
        if EquipInvest is None:
            EquipInvest = self._Mutant__Resources*.2
        self.__ExpSucRate = ErrorRate(EquipInvest, self._Mutant__Resources)
        self._Mutant__Resources -= EquipInvest

            
    def show_BiotechSetting(self):
        '''Report of all properties defined in the biotech experiment.'''
        self.var_Resources = self._Mutant__Resources
        MyVars = [i for i in list(vars(self).keys()) if 'var_' in i]
        for i in range(2): # has to be adjusted to display the Substrate
            print('{}: {}'.format(MyVars[i].replace('var_',''), getattr(self, MyVars[i])))

        
    def show_Library(self):
        '''Report of clones and their performance.'''
        for Clone_ID, Clone_info in self.var_Library.items():
            print("\nClone ID: {}".format(Clone_ID))
            for key in Clone_info:
                print('{}: {}'.format(key, Clone_info[key]))
        
    
    def add_Promoter(self, Clone_ID, Promoter):
        self.var_Library[Clone_ID] = {}
        self.var_Library[Clone_ID]['Promoter_Sequence'] = Promoter
        self.var_Library[Clone_ID]['Promoter_GC-content'] = (Promoter.count('C') + Promoter.count('G')) / len(Promoter)

        
    def Make_MeasurePromoterStrength(self, Clone_ID):
        import random
        ResCost = 100
        if self._Mutant__Resources > ResCost:
            if hasattr(self, 'var_Library'):
                if Clone_ID in self.var_Library:
                    self._Mutant__Resources -= ResCost
                    if random.uniform(0,1) > self._Mutant__ExpSucRate:
                        factor = self._Mutant__InflProStreng          
                        self.var_Library[Clone_ID]['Promoter_Strength'] = round(Help_PromoterStrength(self, Clone_ID) * factor, 2)
                    else:
                        print('Experiment failed, bad equipment.')
                else:
                    print('Error, Clone ID does not exist. Choose existing Clone ID.')
            else:
                print('Error, no promoter library available. Perform a cloning first.')
        else:
            Error_Resources()

                
    def Make_ProductionExperiment(self, Clone_ID, CultTemp, GrowthRate, Biomass, accuracy_Test=.9):
        import numpy as np
        import random
        ResCost = 500
        if self._Mutant__Resources > ResCost: # three resources will be deducted
            # the final experiment can only be performed after at least one sequence has been cloned and tested:
            if hasattr(self, 'var_Library'):
                if Clone_ID in self.var_Library:
                    self._Mutant__Resources -= ResCost
                    if random.uniform(0,1) > self._Mutant__ExpSucRate:
                        # testing whether the determined maximum biomass and the determined maximum growth rate are close to the actual ones
                        if 1 - np.abs(Biomass-self._Mutant__BiomassMax) / self._Mutant__BiomassMax > accuracy_Test and 1 - np.abs(GrowthRate-Help_GrowthConstant(self, self._Mutant__OptTemp)) / Help_GrowthConstant(self, self._Mutant__OptTemp) > accuracy_Test:
                            # Growth rate was only checked, for the calculation the rate resulting from the temperature is used
                            r = Help_GrowthConstant(self, CultTemp)
                            GrowthMax = Growth_Maxrate(self, r, Biomass)
                            self.var_Library[Clone_ID]['Expression_Temperature'] = CultTemp
                            self.var_Library[Clone_ID]['Expression_Biomass'] = Biomass
                            AbsRate = round(GrowthMax * self.var_Library[Clone_ID]['Promoter_Strength'],2)
                            FinalRelRate = round(AbsRate/Calc_MaxExpress(self),2)
                            self.var_Library[Clone_ID]['Expression_Rate'] = FinalRelRate
                            print('{} final vaccine production rate: {}'.format(Clone_ID, FinalRelRate))
                            print('The final production rate is comparable between groups.')
                        else:
                            print('Maximum biomass and/or maximum growth rate are incorrect.')
                    else:
                        print('Experiment failed, bad equipment.')
                else:
                    print('Error, Clone ID does not exist. Choose existing Clone ID.')
            else:
                print('Error, no promoter sequence has been cloned and tested yet. Perform a cloning first and then test the expression with "Make_MeasurePromoterStrength(Clone_ID)".')
                
        else:
            Error_Resources()

    def plot_ReferencePromoterStrength(self):
        '''Function to plot the promoter strength of the optimal sequence additionally as reference.'''
        import matplotlib.pyplot as plt
        
        factor = self._Mutant__InflProStreng
        # Values see init function at the beginning
        if self.var_Host == 'Ecol':
            OptimalPromoterStrength = round(0.057 * factor,2)
        elif self.var_Host == 'Pput':
            OptimalPromoterStrength = round(0.04 * factor,2)
        # plot of maximum Promoter strength together with GC content
        # GC-content is the same for of both optimal sequences.
        plt.plot(0.575, OptimalPromoterStrength, marker = '*', color = 'green', markersize = 10)
            
            
    def Make_TempGrowthExp(self, CultTemps, n=1):
        '''Experiment to determine optimal growth rate. The experiment runs until the maximum biomass is reached.'''
        import numpy as np
        import pandas as pd
        import matplotlib.pyplot as plt
        import time
        import random

        ResCost = 100
        Exp_Duration = 48
        if self._Mutant__Resources > ResCost:
            
            capacity = self._Mutant__BiomassMax
            # the time of the half maximum population (inflection point) is calculated according to here:
            # https://opentextbc.ca/calculusv2openstax/chapter/the-logistic-equation/
            d_mult = 2 # we multiply the inflection point with 'd_mult' to increase cultivation time
            P0 = 0.1
            
            # determine time vector with maximum length:
            OptTemp = self._Mutant__OptTemp
            # Selecting the temperature that is most distant from the optimal temperature
            Temp_tmax = CultTemps[np.argmax(np.absolute(CultTemps-OptTemp))]
            # using the worst temperature to calculate lowest growth rate
            r_tmax = Help_GrowthConstant(self, Temp_tmax)
            # using the worst temperature growth rate to compute longest simulation time, maximum set to 72 h
            duration_tmax = d_mult * 1/r_tmax * np.log((capacity - P0)/P0) + 1
            t_max = np.arange(np.minimum(Exp_Duration, duration_tmax))
            
            # create an empty DataFrame with t_max as first column
            col = []
            col.append('time [h]')
#             col.append('[g/L]:')
            for i in range(len(CultTemps)):
                col.append('exp.{} biomass conc. at {} °C'.format(i, (CultTemps[i])))

            df = pd.DataFrame(np.empty(shape=(len(t_max), len(CultTemps)+1), dtype=float), columns = col)
            df[:len(t_max)] = np.nan 
            new_df = pd.DataFrame({'time [h]': t_max})
            df.update(new_df)
            
            #computing of biomass data and updating of DataFrame
            for i in range(len(CultTemps)):
                if self._Mutant__Resources > ResCost:
                    wait = 0.01 # has to be adjusted, waiting time for loading bar
                    
                    if random.uniform(0,1) > self._Mutant__ExpSucRate: # experiment failure depending on investment to equipment
                        r = Help_GrowthConstant(self, CultTemps[i])
                        # the result can reach very small values, which poses downstream problems, hence the lowest value is set to 0.05
                        if r > 0.05: # under process conditions it might be realistic, source : https://www.thieme-connect.de/products/ebooks/pdf/10.1055/b-0034-10021.pdf
                            duration = Exp_Duration #d_mult * 1/r * np.log((capacity - P0)/P0) + 1
                        else:
                            duration = Exp_Duration
                        t = np.arange(np.minimum(Exp_Duration, duration))
                        
                        # biomass data is calculated according to https://en.wikipedia.org/wiki/Logistic_function
                        mu = capacity / (1 + (capacity-P0) / P0 * np.exp(-r * t))
                        sigma = 0.1*mu
                      
                        exp_TempGrowthExp = [random.normalvariate(mu[k], sigma[k]) for k in range(len(mu))]
                        
                        loading_time = wait * len(t)
                        exp = ' of exp.{} at {} °C'.format(i, (CultTemps[i]))
                        Help_Progressbar(45, loading_time, exp)
                        
                    else:
                        mu = P0
                        sigma = 0.08*mu
                        exp_TempGrowthExp = [random.normalvariate(mu, sigma) for i in range(Exp_Duration)] # if cells haven't grown, the measurement is only continued for 6h
                        
                        loading_time = wait * 7
                        exp = ' of exp.{} at {} °C'.format(i, (CultTemps[i]))
                        Help_Progressbar(45, loading_time, exp)
                
                else:
                    Error_Resources()
                    return
        
                new_df = pd.DataFrame({'exp.{} biomass conc. at {} °C'.format(i, (CultTemps[i])): exp_TempGrowthExp})
                df.update(new_df)
            
                self._Mutant__Resources -= ResCost
            
#             excel_writer = pd.ExcelWriter('Strain_characterization_{}.xlsx'.format(n)) # Export DataFrame to Excel
#             df.to_excel(excel_writer, sheet_name='different temp')
#             excel_writer.close()
            df.to_csv('Strain_characterization_{}.csv'.format(n), index=False)
            
        else:
            Error_Resources()
            return
                
    
    def Make_Cloning(self, Clone_ID, Promoter, Primer, Tm):
        '''Experiment to clone selected promoter. It is displayed whether the experiment was successfull.'''
        import numpy as np
        import random
        
#         if Sequence_ReferenceDistance(Promoter) > .4:
#             return print('Promoter sequence deviates too much from the given structure.')

        ResCost = 200
        if self._Mutant__Resources > ResCost:
            self._Mutant__Resources -= ResCost
            if random.uniform(0,1) < self._Mutant__ExpSucRate:
                print('Experiment failed, bad equipment.')
                return
            
            NaConc = 0.1 # 100 mM source: https://www.genelink.com/Literature/ps/R26-6400-MW.pdf (previous 50 mM: https://academic.oup.com/nar/article/18/21/6409/2388653)
            OptLen = self._Mutant__OptPrLen
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
            
            error_1 = random.uniform(-1,1)*0.1*Primer_Tm_1
            error_2 = random.uniform(-1,1)*0.1*Primer_Tm_2
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
                self.add_Promoter(Clone_ID, Promoter)
            else:
                print('Cloning failed')
                        
        else:
            print('Not enough resources available.')
    
    
    def Choose_Substrate(self, Substrate):
        '''Function to define the C-source for the experiments/predictions'''
        self.var_Substrate = Substrate
            
   
    def ExportExperiments(self):
        '''
        Export of all experiments to a csv file
        '''
        Result_df = NestDict2df(self.var_Library)
        FileName = 'Production_Experiments.csv'
        Result_df.to_csv(FileName, index=True)
        
def Help_PromoterStrength(Mutant, Clone_ID, Predict_File=None, Similarity_Thresh=.4):
    '''Expression of the recombinant protein.
        Arguments:
            Mutant: class, contains optimal growth temperature, production phase
            Clone_ID: Clone with defined promoter for which express
            Predict_File: string, address of regression file
        Output: 
            Expression: float, expression rate
    '''
    import os
    import numpy as np
    import joblib
    import pickle
    
    if Sequence_ReferenceDistance(Mutant.var_Library[Clone_ID]['Promoter_Sequence']) > Similarity_Thresh:
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

        X_Test = np.array(list_onehot(np.delete(list_integer(Mutant.var_Library[Clone_ID]['Promoter_Sequence']),Positions_removed, axis=0))).reshape(1,-1)  
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
    # calculation of the growth rate constant, by picking the activity from a normal distribution
    
    growth_rate_const = r_pdf.pdf(CultTemp) / r_pdf.pdf(OptTemp)
    
    return growth_rate_const
    

def Growth_Maxrate(Mutant, growth_rate_const, Biomass):
    '''The function calculates the maximum slope during growth.
        Arguments:
            Mutant: class, contains maximum biomass concentration as carrying capacity
            growth_rate_const: float, maximum growth rate constant
        Output:
            growth_rate_max: float, maximum growth rate
    '''
    # biomass checks
#     if Biomass > Mutant._Mutant__BiomassMax or not Biomass:
#         print('Error, no biomass was set or unexpected value or the maximum possible biomass was exceeded. Enter a value for the biomass again.')        
    
    # Equation for calculating the maximum slope
    # https://www.tjmahr.com/anatomy-of-a-logistic-growth-curve/
    GrowthMax = Biomass * growth_rate_const / 4
    
    return GrowthMax
    
    

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
    
    
def Help_SwitchComplementary(argument):
    switcher = {
        'T': 'A',
        'A': 'T',
        'C': 'G',
        'G': 'C'
    }
    return switcher.get(argument)

def Error_Resources():
    print('Not enough resources available.')
    
    
def Plot_ExpressionRate():
    '''function to plot the expression rate as a function of growth rate and promoter strength. The aim was to find out how both influencing variables are in the same order of magnitude.     
    assumption: The values of the promoter strength (y) are in range from 0.001 to 0.025 in case of *P. putida* respectively in range from 0.001 to 0.05 in case of E. coli.    
    The values of the growth rate (x) are in range from 0 to 1 on the basis of the standardisation.    
    The factor was determined using the values for P. putida. Accordingly, for E. coli the values for promoter strength and expression rate should be twice as high at the end.'''
    import numpy as np
    import matplotlib.pyplot as plt

    x = np.linspace(0, 1, 50)

    ymin = 30*np.linspace(0.001, 0.025, 50)
    ym = 40*np.linspace(0.001, 0.025, 50) # factor = 1:0,025, so that x and y are in the same order of magnitude
    ymax = 50*np.linspace(0.001, 0.025, 50)

    X, Ymin = np.meshgrid(x, ymin)
    X, Ym = np.meshgrid(x, ym)
    X, Ymax = np.meshgrid(x, ymax)
    
    # 3D plot to to find out the connection and visualize the influence of the factors
    fig = plt.figure(figsize = (7,6), dpi = 120)
    ax = plt.axes(projection='3d')

    Z = np.multiply(X,Ymin)
    ax.contour3D(X, Ymin, Z, 20, cmap='binary')
    Z = np.multiply(X,Ym)
    ax.contour3D(X, Ym, Z, 20)
    Z = np.multiply(X,Ymax)
    ax.contour3D(X, Ymax, Z, 20, cmap='inferno')

    ax.set_xlabel('normalized growth rate [-]')
    ax.set_ylabel('promoter strength [-]')
    ax.set_zlabel('expression rate [-]');
    
    '''As can be seen in the plot, the factor by which the promoter strength is multiplied does not change
        the influence of this strength.
        The factor can be used to influence the range of the promoter strength and thus the range
        of the final expression rate.'''
    
def ErrorRate(Invest, ResTotal= 5000, relKM=.005, Vmax=.9):
    '''
    Calculate the experimental error rate based on the investment to equipment. Based on inverse Michaelis-Menten equation mit max error 0.9 and min error 0.1.
    '''
    KM = relKM*ResTotal
    myMax = 1
    return myMax - (Vmax * Invest) / (KM + Invest) 


def Calc_MaxExpress(self):
    '''Function to calculate the maximum possible expression rate.'''
    BiomassMax = self._Mutant__BiomassMax
    OptTemp = self._Mutant__OptTemp
    factor = self._Mutant__InflProStreng
    # Values see init function at the beginning
    if self.var_Host == 'Ecol':
        MaximumPromoterStrength = round(0.057 * factor,2)
    elif self.var_Host == 'Pput':
        MaximumPromoterStrength = round(0.04 * factor,2)
    r = Help_GrowthConstant(self, OptTemp)
    GrowthMax = Growth_Maxrate(self, r, BiomassMax)
    return round(GrowthMax * MaximumPromoterStrength,2)


def NestDict2df(mydict):
    '''
    Converts a nested dictionary into a dataframe
    '''
    import pandas as pd
    
    return pd.concat([pd.DataFrame.from_dict(mydict[who], orient='index', columns=[who]) for who in mydict.keys()], axis=1).transpose()