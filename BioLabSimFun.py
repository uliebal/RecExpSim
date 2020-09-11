# Biotechnology Laboratory Simulator Functions
# print(vars(myhost))
 
class Host:
    '''The 'Host' class stores all information about the organism and the integrated recombinant protein.''' 
    
#     from BioLabSimFun import SequenceRandomizer_Single
    from random import randint
    # random assignment of the production phase, either during growth phase or stationary phase
    __ProdPhase = 'exponential' if randint(0,1)==0 else 'stationary'
    # resources, e.g. money, for conducting tests
    __Resources = 40
    __BiomassMax = None
    
    def __init__(self, Host, MetNet = False):
        import os
        from random import randint
        self.var_Host = Host
        self.var_Resources = self._Host__Resources
        self.var_Substrate = None
        # Library variable containing details to the different tested mutants
        self.var_Library = {}
        # Strain collection containing manipulated genomes
        self.info_Strain = {}
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
        # initiating metabolic network

    
    def show_BiotechSetting(self):
        '''Report of all properties defined in the biotech experiment.'''
        self.var_Resources = self._Host__Resources
        MyVars = [i for i in list(vars(self).keys()) if 'var_' in i]
        for i in range(len(MyVars)): # has to be adjusted to display the Substrate
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
        if self._Host__Resources > 0:
            if hasattr(self, 'var_Library'):
                if Clone_ID in self.var_Library:
                    factor = self._Host__InflProStreng 
                    Sequence = Host.var_Library[Clone_ID]['Promoter_Sequence']
                    self.var_Library[Clone_ID]['Promoter_Strength'] = round(Help_PromoterStrength(self.var_Host, Sequence) * factor, 2)
                    self._Host__Resources -= 1
                else:
                    print('Error, Clone ID does not exist. Choose existing Clone ID.')
            else:
                print('Error, no promoter library available. Perform a cloning first.')
        else:
            Error_Resources()

                
    def Make_ProductionExperiment(self, Clone_ID, CultTemp, GrowthRate, Biomass, accuracy_Test=.9):
        import numpy as np
        if self._Host__Resources > 2: # three resources will be deducted
            # the final experiment can only be performed after at least one sequence has been cloned and tested:
            if hasattr(self, 'var_Library'):
                if Clone_ID in self.var_Library:
                    # testing whether the determined maximum biomass and the determined maximum growth rate are close to the actual ones
                    if 1 - np.abs(Biomass-self._Host__BiomassMax) / self._Host__BiomassMax > accuracy_Test and 1 - np.abs(GrowthRate-Help_GrowthConstant(self, self._Host__OptTemp)) / Help_GrowthConstant(self, self._Host__OptTemp) > accuracy_Test:
                        # Growth rate was only checked, for the calculation the rate resulting from the temperature is used
                        r = Help_GrowthConstant(self, CultTemp)
                        GrowthMax = Growth_Maxrate(self, r, Biomass)
                        self.var_Library[Clone_ID]['Expression_Temperature'] = CultTemp
                        self.var_Library[Clone_ID]['Expression_Biomass'] = Biomass
                        self.var_Library[Clone_ID]['Expression_Rate'] = round(GrowthMax * self.var_Library[Clone_ID]['Promoter_Strength'],2)
                        self._Host__Resources -= 3
                    else:
                        print('Maximum biomass and/or maximum growth rate are incorrect.')
                else:
                    print('Error, Clone ID does not exist. Choose existing Clone ID.')
            else:
                print('Error, no promoter sequence has been cloned and tested yet. Perform a cloning first and then test the expression with "Make_MeasurePromoterStrength(Clone_ID)".')
                
        else:
            Error_Resources()


    def show_TargetExpressionRate(self):
        '''Function to calculate the maximum possible expression rate and to tell the students what the minimum rate should be.'''
        BiomassMax = self._Host__BiomassMax
        OptTemp = self._Host__OptTemp
        factor = self._Host__InflProStreng
        # Values see init function at the beginning
        if self.var_Host == 'Ecol':
            MaximumPromoterStrength = round(0.057 * factor,2)
        elif self.var_Host == 'Pput':
            MaximumPromoterStrength = round(0.04 * factor,2)
        r = Help_GrowthConstant(self, OptTemp)
        GrowthMax = Growth_Maxrate(self, r, BiomassMax)
        achievExpRate = round(0.75*GrowthMax * MaximumPromoterStrength,2)
        print('At least an expression rate of {} should be achieved by the production experiment.'.format(achievExpRate))
        
        
    def plot_ReferencePromoterStrength(self):
        '''Function to plot the promoter strength of the optimal sequence additionally as reference.'''
        import matplotlib.pyplot as plt
        
        factor = self._Host__InflProStreng
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

        if self._Host__Resources > 0:
            
            capacity = self._Host__BiomassMax
            # the time of the half maximum population (inflection point) is calculated according to here:
            # https://opentextbc.ca/calculusv2openstax/chapter/the-logistic-equation/
            d_mult = 2 # we multiply the inflection point with 'd_mult' to increase cultivation time
            P0 = 0.1
            
            # determine time vector with maximum length:
            OptTemp = self._Host__OptTemp
            # Selecting the temperature that is most distant from the optimal temperature
            Temp_tmax = CultTemps[np.argmax(np.absolute(CultTemps-OptTemp))]
            # using the worst temperature to calculate lowest growth rate
            r_tmax = Help_GrowthConstant(self, Temp_tmax)
            # using the worst temperature growth rate to compute longest simulation time, maximum set to 72 h
            duration_tmax = d_mult * 1/r_tmax * np.log((capacity - P0)/P0) + 1
            t_max = np.arange(np.minimum(73, duration_tmax))
            
            # create an empty DataFrame with t_max as first column
            col = []
            col.append('time [h]')
            col.append('[g/L]:')
            for i in range(len(CultTemps)):
                col.append('exp.{} biomass conc. at {} °C'.format(i+1, (CultTemps[i])))

            df = pd.DataFrame(np.empty(shape=(len(t_max), len(CultTemps)+2), dtype=float), columns = col)
            df[:len(t_max)] = np.nan 
            new_df = pd.DataFrame({'time [h]': t_max})
            df.update(new_df)
            
            #computing of biomass data and updating of DataFrame
            for i in range(len(CultTemps)):
                if self._Host__Resources > 0:
                    wait = 0.01 # has to be adjusted, waiting time for loading bar
                    
                    if random.uniform(0,1) > 0.1: # in 10% of cases the culture does not grow (failure of the experiment)
                        r = Help_GrowthConstant(self, CultTemps[i])
                        # the result can reach very small values, which poses downstream problems, hence the lowest value is set to 0.05
                        if r > 0.05: # under process conditions it might be realistic, source : https://www.thieme-connect.de/products/ebooks/pdf/10.1055/b-0034-10021.pdf
                            duration = d_mult * 1/r * np.log((capacity - P0)/P0) + 1
                        else:
                            duration = 7
                        t = np.arange(np.minimum(73, duration))
                        
                        # biomass data is calculated according to https://en.wikipedia.org/wiki/Logistic_function
                        mu = capacity / (1 + (capacity-P0) / P0 * np.exp(-r * t))
                        sigma = 0.1*mu
                      
                        exp_TempGrowthExp = [random.normalvariate(mu[k], sigma[k]) for k in range(len(mu))]
                        
                        loading_time = wait * len(t)
                        exp = ' of exp.{} at {} °C'.format(i+1, (CultTemps[i]))
                        Help_Progressbar(45, loading_time, exp)
                        
                    else:
                        mu = P0
                        sigma = 0.08*mu
                        exp_TempGrowthExp = [random.normalvariate(mu, sigma) for i in range(7)] # if cells haven't grown, the measurement is only continued for 6h
                        
                        loading_time = wait * 7
                        exp = ' of exp.{} at {} °C'.format(i+1, (CultTemps[i]))
                        Help_Progressbar(45, loading_time, exp)
                
                else:
                    Error_Resources()
                    return
        
                new_df = pd.DataFrame({'exp.{} biomass conc. at {} °C'.format(i+1, (CultTemps[i])): exp_TempGrowthExp})
                df.update(new_df)
            
                self._Host__Resources -= 1
            
            excel_writer = pd.ExcelWriter('Strain_characterization_{}.xlsx'.format(n)) # Export DataFrame to Excel
            df.to_excel(excel_writer, sheet_name='different temp')
            excel_writer.close()
            
        else:
            Error_Resources()
            return
                
    
    def Make_Cloning(self, Clone_ID, Promoter, Primer, Tm):
        '''Experiment to clone selected promoter. It is displayed whether the experiment was successfull.'''
        import numpy as np
        import random
        
        if Sequence_ReferenceDistance(Promoter) > .4:
            return print('Promoter sequence deviates too much from the given structure.')
        
        if self._Host__Resources > 0:
            
            NaConc = 0.1 # 100 mM source: https://www.genelink.com/Literature/ps/R26-6400-MW.pdf (previous 50 mM: https://academic.oup.com/nar/article/18/21/6409/2388653)
            OptLen = self._Host__OptPrLen
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
                
            self._Host__Resources -= 1
        
        else:
            print('Not enough resources available.')
    
    
#     def Choose_Substrate(self, Substrate):
#         '''Function to define the C-source for the experiments/predictions'''
#         self.var_Substrate = Substrate
            
class Strain:
    '''
    The 'strain' class stores information of genome and other associated metabolic information
    '''
    def __init__(self, Host = 'Ecol', GenomeSize = 500, GCcont = .6):
        import os
        print('Initiating metabolic network')
        ModelPath = os.path.join('Models','e_coli_core.xml')
        self.var_Model = Help_LoadCobra(Path = ModelPath)
        self.info_GenesDF = Help_GeneAnnotator(Host, self.var_Model) # self.__GenesDF = 
        self.info_GenesDF['Expr2Flux'] = Help_Expr2Flux(self.info_GenesDF) # self.__GenesDF = 
        self.info_Genome = Help_GenomeGenerator(self.info_GenesDF, GenomeSize, GCcont)
        self.var_GenomeSize = len(self.info_Genome)
        self.var_GCcont = round((self.info_Genome.count('G') + self.info_Genome.count('C'))/self.var_GenomeSize,2)
#             else:
#                 print('Invalid organism for metabolic simulation.')
        
def Help_PromoterStrength(Host, Sequence, Scaler=100, Similarity_Thresh=.4, Predict_File=None):
    '''Expression of the recombinant protein.
        Arguments:
            Host:       class, contains optimal growth temperature, production phase
            Sequence:     string, Sequence for which to determine promoter strength
            Scaler:       int, multiplied to the regression result for higher values
            Predict_File: string, address of regression file
        Output: 
            Expression: float, expression rate
    '''
    import os
    import numpy as np
    import joblib
    import pickle
    
    if Sequence_ReferenceDistance(Sequence) > Similarity_Thresh:
        Expression = 0
    else:
        if Predict_File!=None:
            Regressor_File = Predict_File
        else:    
            Data_Folder = 'ExpressionPredictor'
            if Host == 'Ecol':
                Regressor_File = os.path.join(Data_Folder,'Ecol-Promoter-predictor.pkl')
                Add_Params = os.path.join(Data_Folder,'Ecol-Promoter-AddParams.pkl')
#                 Scaler_DictName = 'Ecol Promoter Activity_Scaler'
            elif Host == 'Pput':
                Regressor_File = os.path.join(Data_Folder,'Ptai-Promoter-predictor.pkl')
                Add_Params = os.path.join(Data_Folder,'Ptai-Promoter-AddParams.pkl')
#                 Scaler_DictName = 'Ptai Promoter Activity_Scaler'
            else:
                print('Non-recognized host name. Rename host to either "Ecol" or "Pput."')

        Predictor = joblib.load(Regressor_File)
        Params = pickle.load(open(Add_Params, 'rb'))
        Positions_removed = Params['Positions_removed']
#         Expr_Scaler = Params[Scaler_DictName]

        X = np.array(list_onehot(np.delete(list_integer(Sequence),Positions_removed, axis=0))).reshape(1,-1)  
        GC_cont = (Sequence.count('G') + Sequence.count('C'))/len(Sequence)
        X = np.array([np.append(X,GC_cont)])
        Y = Predictor.predict(X)
        Expression = round(float(Y)*Scaler,3)

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

def Help_GrowthConstant(Host, CultTemp, var=5):
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
    
    OptTemp = Host._Host__OptTemp
    r_pdf = norm(OptTemp, var)
    # calculation of the growth rate constant, by picking the activity from a normal distribution
    
    growth_rate_const = r_pdf.pdf(CultTemp) / r_pdf.pdf(OptTemp)
    
    return growth_rate_const
    

def Growth_Maxrate(Host, growth_rate_const, Biomass):
    '''The function calculates the maximum slope during growth.
        Arguments:
            Host: class, contains maximum biomass concentration as carrying capacity
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


def Help_ExportToExcel(Host, FileName, sheet_name, x_values, y_values,
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
    
def Help_LoadCobra(Path=False):
    '''
    This function loads cobra models from input path, or by default the Ecoli core model. Either from a local copy or dowloading it from the Bigg database.
    
    Input:
        Path (optional): string, relative path address to local copy
        
    Output:
        Model:           cobraby model, E. coli core reactions
    '''
    
    import os
    from cobra.io import read_sbml_model    
    
    if Path:
        Model = read_sbml_model(Path)
    else:
        os.system('wget http://bigg.ucsd.edu/static/models/e_coli_core.xml')
        # generating cobra variable from SBML/xml/mat file
        Path = 'e_coli_core.xml'
        Model = read_sbml_model(Path)
        
    return Model


def Help_GeneAnnotator(Host, Model):
    '''
    Function to characterize genes.
    Input
        Model:    cobra model
    Output
        Genes_df: dataframe, gene name from enzyme id in model, expression strength, promoter sequence, ORF, flux
    '''
    import pandas as pd
    import multiprocessing
    from joblib import Parallel, delayed
    
    num_cores = multiprocessing.cpu_count()
#     use_core = min([num_cores, n])
  
    Result = Parallel(n_jobs=num_cores)(delayed(make_GeneJoiner)(Host, Model, myRct.id) for myRct in Model.reactions)
    Genes_df = pd.DataFrame(Result)
    # adding flux values
    # setup of flux boundaries. For the reference boundary changes are set to 'False', 
    # for mutant strains, ractions with altered promoter sequence will change enzyme levels and boundaries must be changed accordingly, their variable is 'True'
#     myModel = make_AdaptModel(Model, Strain)
    Fluxes = Model.optimize()
    Genes_df['Fluxes'] = Fluxes.fluxes.values
    
    return Genes_df

def Help_GenomeGenerator(GenesDF, GenomeSize, GCcont):
    '''
    Constructs whole genome with interspersed genes.
    '''
    import numpy as np
    import random
    # combining promoter and ORF
    Genes = [''.join([myProm, myORF]) for myProm, myORF in zip(GenesDF['Promoter'].values,GenesDF['ORF'].values)]
    Genes_List = [Convert(Gene) for Gene in Genes]
    # generating background genome sequence
    Genome_Bckgd = make_GenomeBckgd(GenomeSize, GCcont)
    # determining position for gene insertion
    Gene_Positions = np.sort(random.sample(range(len(Genome_Bckgd)),len(Genes)))
    # breaking the background genome in nested lists at gene positions
    Genome_Tmp = make_NestedList(Genome_Bckgd, Gene_Positions)
    # Now inserting the genes
    Gtmp = [np.concatenate([mbg,bed]) for mbg,bed in zip(Genome_Tmp[:-1],Genes_List)]
    Gtmp = np.concatenate([Gtmp,Genome_Tmp[-1]])
    Genome = ''.join([''.join(elm) for elm in Gtmp])    
    
    return Genome

def Help_StrainCharacterizer(Host, GenesDF, Genome_WT, Genome_MT, Model):
    '''
    Function to scan a manipulated genome for changes in gene expression of enzymes.
    Input:
        Mutant:         parent class
        StrainID
    Output:
        RctNew_df:      dataframe, containing reactions with updated expression
    '''
    import pandas as pd
    import multiprocessing
    from joblib import Parallel, delayed
    
    num_cores = multiprocessing.cpu_count()
#     use_core = min([num_cores, n])
  
    Result = Parallel(n_jobs=num_cores)(delayed(make_UpdateExpression)(Host, GenesDF, Genome_WT, Genome_MT, myRct.id) for myRct in Model.reactions)
    RctNew_df = pd.DataFrame(Result)

    return RctNew_df
    
def Help_Expr2Flux(GenesDF):
    '''
    Function to correlate expression strength with flux.
    Input:
        GenesDF:        DataFrame, details of enzymes
    Output:
        Corr_ExprFlux:  array, correlation factor
    '''
    Corr_ExprFlux = GenesDF['Fluxes'].values/GenesDF['Expression'].values
    
    return Corr_ExprFlux

def make_GeneJoiner(Host, Model, RctID):
    '''
    Determines promoter activity and combines it with enzyme id and ORF.
    Output
        Gene_Info:    dictionary, gene id, gene expression, promoter, ORF
    '''

    Gene_ORF = make_ORF(Model)
    Gene_Promoter = make_Promoter()
    Gene_Activity = Help_PromoterStrength(Host, Gene_Promoter, Similarity_Thresh=.8)
    
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
    import random

    SeqNestList = [random.choices([Letter for Nest in random.choices([['G','C'],['A','T']], weights=[GCcont, 1-GCcont]) for Letter in Nest]) for _x in range(GenomeSize)]
    Genome_Bckgd = ''.join(Letter for Nest in SeqNestList for Letter in Nest)

    return Genome_Bckgd

def make_ORF(Model):
    '''
    Function to generate a single ORF. Triplet frequencies are taken from E.coli:
    https://openwetware.org/wiki/Escherichia_coli/Codon_usage
    ORF starts always with 'ATG'.
    Input
        CodonTriplets:     dataframe, automatic load, base triplets, ID (e.g. 'Stop', 'Met'), frequency in percent
        Mutant:            class, subfield var_Model.reactions is used to count the enzyme number to derive minimum sequence length
    Output
        Gene_ORF:          string, open reading frame of enzyme
    '''
    import random
    import numpy as np
    import pandas as pd

    # coding sequence construction
    # first we determine the minimum coding gene length of nucleotides to distinguish the enzymes in the model
    Enzyme_Number = len(Model.reactions)
    Gen_Minimum = np.ceil(np.log2(Enzyme_Number))

    # we want to represent codon triplicates, we calculate the next highest divisor of three
    Gene_Length = int(np.ceil(Gen_Minimum/3))
    CodonFile = 'CodonTriplets.csv'
    CodonTriplets = pd.read_csv(CodonFile, delimiter=';', skipinitialspace=True)
    CodonStop = CodonTriplets[['Stop' in s for s in CodonTriplets['Name']]].reset_index()
    CodonCoding = CodonTriplets.drop(CodonStop.index).reset_index()
    Gene_ORF = [random.choices(CodonCoding['Triplet'], weights=CodonCoding['Percent']) for CodonId in range(Gene_Length)]
    Gene_Stop = random.choices(CodonStop['Triplet'], weights=CodonStop['Percent'])
    Gene_ORF.append(Gene_Stop)
    Gene_ORF.insert(0,['ATG'])
    Gene_ORF = ''.join([Letter for Nest in Gene_ORF for Letter in Nest])

    return Gene_ORF

def make_Promoter(RefFile=False, WeightFile=False):
    '''
    Generating a single promoter. The promoter sequences are selected from the exploration space used to train the regressor for promoter strength. The algorithm needs to know which nucleotides at each position were sampled strongly enough for predictions. Positions not sufficiently well sampled are filled with entries from a reference sequence.
    Input
        RefFile:               string, path to reference sequence in txt-file, used to fill positions with weak sampling
        WeightFile:            string, path to pickle with dataframe, columns represent bases (A-C-G-T), rows represent positions, index is position relative to start codon (e.g. -40, -1)
    Output
        Gene_Promoter:        string, sequence of the promoter, same length as reference sequence
    '''
    import numpy as np
    import pickle
    import os
    import random
    # generating promoter sequence

    # the reference sequence contains the most common tested nucleotides at each position
    if not RefFile:   
        RefFile = os.path.join('Models','RefSeq.txt')
    with open(RefFile) as f:
        RefSeq = f.read()
    
    # loading information that determines the exploratory space of the regressor
    # weight file has boolean representations for each base on each position whether it was part of the training set or not
    if not WeightFile:
        WeightFile = os.path.join('Models','NucleotideWeightTable.pkl')
    with open(WeightFile, 'rb') as handle:
        Nucleotides_Weight = pickle.load(handle)
    Nucleotides_Weight.index = Nucleotides_Weight.index+40

    # setup of default sequence with reference
    SeqDef = np.asarray([Letter for Letter in RefSeq])

    # changing positions within the exploratory space of the regressor
    Bases = ['A','C','G','T']
    for idx in range(len(Nucleotides_Weight)):
        np.put(SeqDef,Nucleotides_Weight.index[idx],random.choices(Bases,Nucleotides_Weight.iloc[idx].values))
    Gene_Promoter = ''.join(SeqDef)
    
    return Gene_Promoter

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

def make_UpdateExpression(Host, GenesDF, Genome_WT, Genome_Mut, RctID):
    '''
    Function to evaluate if expression for the input gene is different compared to the reference strain. 
    Input:
        Mutant:         parent class, subfield _Mutant_GenesDF and info_Genome used.
        StrainID:       string, identifier for expressed genome in subclass var_strains
        RctID:          string, identifier for the reaction to be tested
    Output:
        OutDict:        dictionary containing 'RctFlag' boolean, if True, the expression relative to the reference has changed; 'Gene_Activity' float, expression strength; 'RctID' string, identifier for the reaction tested
    '''
    
    import numpy as np
    
    # Extraction of ORF and reference promoter sequence of input reaction from reference annotation database
    Rct_Idx = GenesDF[GenesDF['RctID']==RctID].index.values
    Rct_ORF = ''.join(str(n) for n in GenesDF['ORF'].iloc[Rct_Idx].to_list())
    Rct_RefProm = GenesDF['Promoter'].iloc[Rct_Idx].values

    # Extraction of reaction association promoter sequence from input genome
    Gene_ORF = Genome_WT.find(Rct_ORF)
    Gene_Promoter = Genome_Mut[Gene_ORF-40:Gene_ORF:1]

    # Comparison of current promoter with reference promoter
    # If promoters differ, then new expression strength is calculated
    if Gene_Promoter != Rct_RefProm:
#         print('New promoter.')
        RctFlag = True
        Gene_Activity = float(Help_PromoterStrength(Host, Gene_Promoter, Similarity_Thresh=.8))
    else:
        RctFlag = False
        Gene_Activity = float(GenesDF['Expression'].iloc[Rct_Idx].values)
#     print(Gene_Activity)
    OutDict = {'RctFlag':RctFlag, 'Activity':Gene_Activity, 'RctID':RctID}
    
    return OutDict

def make_AdaptModel(StrainID):
    '''
    Function to adapt the Cobra-model boundaries according to the expression strength of the genome promoters.
    Input:
        :         parent class, subfield _Mutant_GenesDF and info_Genome used.
        StrainID:       string, identifier for expressed genome in subclass var_strains
    Output:
        Model:          cobraby model, E. coli core reactions with updated boundary conditions
    '''
    
    
    return Model