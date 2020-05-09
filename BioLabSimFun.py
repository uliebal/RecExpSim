# Biotechnology Laboratory Simulator Functions

class Mutant:
    '''The 'Mutant' class stores all information about the organism and the integrated recombinant protein.''' 
    
#     from BioLabSimFun import SequenceRandomizer_Single
    from random import randint
    # random assignment of the production phase, either during growth phase or stationary phase
    __ProdPhase = 'exponential' if randint(0,1)==0 else 'stationary'
    # resources, e.g. money, for conducting tests
    __Resources = 50
    __BiomassMax = None
    
    def __init__(self, Host):
        from random import randint
        self.var_Host = Host
#         self.var_Promoter = []
        self.var_Resources = self._Mutant__Resources
        # factor which influences the range of the promoter strength, randomly assigned
        self.__InflProStreng = randint(30,50) # explanation see workflow 
        # optimal growth temperature, randomly assigned
        self.__OptTemp = randint(25,40) # unit: degree celsius, source 1: https://refubium.fu-berlin.de/bitstream/handle/fub188/7617/02_2_1_Literatur.pdf?sequence=3&isAllowed=y, source 2:https://www.baua.de/DE/Angebote/Rechtstexte-und-Technische-Regeln/Regelwerk/TRBA/pdf/Pseudomonas-putida.pdf?__blob=publicationFile&v=2
        # optimal Primer length, randomly assigned
        self.__OptPrLen = randint(16,28) # unit: nt, source: https://link.springer.com/article/10.1007/s10529-013-1249-8
        # maximum biomass concentration, can be adjusted later, now randomly set
        if self.var_Host == 'Ecol':
            self.__BiomassMax = randint(10,160) # unit: in gDCW/l, source: https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=2&cad=rja&uact=8&ved=2ahUKEwjzt_aJ9pzpAhWGiqQKHb1jC6MQFjABegQIAhAB&url=https%3A%2F%2Fwww.repo.uni-hannover.de%2Fbitstream%2Fhandle%2F123456789%2F3512%2FDissertation.pdf%3Fsequence%3D1&usg=AOvVaw2XfGH11P9gK2F2B63mY4IM
        elif self.var_Host == 'Pput':
            self.__BiomassMax = randint(30,100) # unit: in gDCW/l, source 1: https://onlinelibrary.wiley.com/doi/pdf/10.1002/bit.25474, source 2: https://link.springer.com/article/10.1385/ABAB:119:1:51
    
    def show_BiotechSetting(self):
        '''Report of all properties defined in the biotech experiment.'''
        self.var_Resources = self._Mutant__Resources
        MyVars = [i for i in list(vars(self).keys()) if 'var_' in i]
        for i in range(len(MyVars)):
            print('{}: {}'.format(MyVars[i].replace('var_',''), getattr(self, MyVars[i])))

    #def show_Primer_DeviationOptimalLength(self, PrLen):
        #import numpy as np
        #OptLen = self._Mutant__OptPrLen
        #Devi = (np.absolute(OptLen - PrLen)/OptLen)*100
        #print('The deviation from the optimum length is {} %.'.format(Devi.round(2)))
        #if PrLen <= 30:
            #print('The length of the primer does not exceed 30 nt.')
        #else:
            #print('The primer is too long (>30 nt).')
    
    def add_Promoter(self, Promoter):
        self.var_Promoter = Promoter
        self.var_GC_content = (self.var_Promoter.count('C') + self.var_Promoter.count('G')) / len(self.var_Promoter)

    def __add_RandomPromoter(self):
        self.var_Promoter = SequenceRandomizer_Single()
        self.var_GCcontent = (self.var_Promoter.count('C') + self.var_Promoter.count('G')) / len(self.var_Promoter)
        
    def Make_MeasurePromoterStrength(self):
        if not hasattr(self, 'var_Promoter'):
            print('Error, no promoter added. Perform a cloning first.')
            return
        if self._Mutant__Resources > 0:
            factor = self._Mutant__InflProStreng          
            self.var_PromoterStrength = round(Help_Expression(self) * factor, 2)
            self._Mutant__Resources -= 1
        else:
            print('Not enough resources available.')
            self.var_PromoterStrength = None
        
    def Make_ProductionExperiment(self):
        # so that the final experiment can only be performed after at least one sequence has been cloned and tested:
        if not hasattr(self, 'var_PromoterStrength'):
            print('Error, no promoter sequence has been cloned and tested yet. Perform a cloning first and then test the expression with "Make_MeasurePromoterStrength()".')
            return
        
        if self._Mutant__Resources > 1: # two resources will be deducted
            # first the selected promoter sequence must be entered so that the experiment can be performed
            
            ReEntry = 1
            while ReEntry:
                ReEntry = 0
                try:
                    Promoter = input('Choose a promoter sequence: ')
                    if not Promoter:
                        raise ValueError('empty string')
                except ValueError:
                    print('Error, none of the produced recombinant strains with the desired promoter was selected. Choose a promoter sequence first.')
                    ReEntry = 1
            
            self.add_Promoter(Promoter)
            # first, the required promoter strength is measured again
            self.Make_MeasurePromoterStrength()
            self._Mutant__Resources -= 1
            
            # input of the temperature
            # if no tempertaure is set, the optimal one is used
            try:
                CultTemp = int(input('temperature [°C] for the experiment: '))
            except ValueError:
                CultTemp = self._Mutant__OptTemp
                print('No temperature was set or unexpected value, therefore the optimal temperature was used.')
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
        import random
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
                #random.seed()
                if self._Mutant__Resources > 0:
                    wait = 0.01 # has to be adjusted, waiting time for loading bar
                    
                    if random.uniform(0,1) > 0.1: # in 10% of cases the culture does not grow (failure of the experiment)
                        #print('Temp Test {}:', CultTemps[i])
                        r = Help_GrowthConstant(self, CultTemps[i])
                        # the result can reach very small values, which poses downstream problems, hence the lowest value is set to 0.05
                        if r > 0.05: # under process conditions it might be realistic, source : https://www.thieme-connect.de/products/ebooks/pdf/10.1055/b-0034-10021.pdf
                            duration = d_mult * 1/r * np.log((capacity - P0)/P0) + 1
                        else:
                            duration = 7
                        t = np.arange(np.minimum(73, duration))
                    
                        mu = capacity / (1 + (capacity-P0) / P0 * np.exp(-r * t))
                        sigma = 0.1*mu
                      
                        exp_TempGrowthExp = random.normalvariate(mu, sigma)
                        
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
                    print('Not enough resources available.')
                    return
        
                new_df = pd.DataFrame({'exp.{} biomass conc. at {} °C'.format(i+1, (CultTemps[i])): exp_TempGrowthExp})
                df.update(new_df)
            
                self._Mutant__Resources -= 1
            
            excel_writer = pd.ExcelWriter('Strain_characterization.xlsx') # Export DataFrame to Excel
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
        
    
    def Make_Cloning(self, Primer, Tm, Promoter):
        '''Experiment to clone selected promoter. It is displayed whether the experiment was successfull.'''
        import numpy as np
        import random
        
        if self._Mutant__Resources > 0:
            
            NaConc = 5e-02 # 50 mM source: https://academic.oup.com/nar/article/18/21/6409/2388653
            OptLen = self._Mutant__OptPrLen
            AllowDevi = 0.2 # allowed deviation
            Primer_Length = len(Primer)
            Primer_nC = Primer.count('C')
            Primer_nG = Primer.count('G')
            Primer_GC_content = (Primer_nC + Primer_nG) / Primer_Length
            
            Primer_Tm = 81.5 + 16.6*np.log10(NaConc) + 0.41*Primer_GC_content - 675/Primer_Length # source: https://core.ac.uk/download/pdf/35391868.pdf#page=190
            # Product_Tm = 0.41*(Primer_GC_content) + 16.6*np.log10(NaConc) - 675/Product_Length
            # Ta_Opt = 0.3*Primer_Tm + 0.7*Product_Tm - 14.9
            # source Product_Tm und Ta: https://academic.oup.com/nar/article/18/21/6409/2388653
            # Product_Length wäre die Länge des Promoters (40)? zu klein -> negative Zahl kommt raus für Product_Tm
            
            error = random.uniform(-1,1)*0.1*Primer_Tm
            Primer_Tm_err = error + Primer_Tm
            
            DeviLen = np.absolute(OptLen - Primer_Length)/OptLen
            DeviTm = np.absolute(Primer_Tm_err - Tm)/Primer_Tm_err
            
            #create the complementary sequence of the primer to check for mistakes:
            PrimerComp = ""
            for base in Primer:
                PrimerComp = PrimerComp + Help_SwitchComplementary(base)
            
            if DeviLen <= AllowDevi and DeviTm <= AllowDevi and Primer_Length <= 30 and PrimerComp == Promoter[:len(Primer)]:
                print('Cloning was successfull.')
                self.add_Promoter(Promoter)
                #exp_Cloning = (1 - np.absolute(Primer_Tm_err - Tm)/Primer_Tm_err) * 100
                #print(f'The efficiency of cloning is {exp_Cloning.round(2)} %.')
            
            else:
                print('Cloning failed')
                
            self._Mutant__Resources -= 1
        
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
    # calculation of the growth rate constant, by picking the activity from a normal distribution
    
    growth_rate_const = r_pdf.pdf(CultTemp) / r_pdf.pdf(OptTemp)
    
    return growth_rate_const
    

def Growth_Maxrate(Mutant, growth_rate_const):
    '''The function calculates the maximum slope during growth.
        Arguments:
            Mutant: class, contains maximum biomass concentration as carrying capacity
            growth_rate_const: float, maximum growth rate constant
        Output:
            growth_rate_max: float, maximum growth rate
    '''
    # input of the biomass
    # if no biomass is set or if the maximum is exceeded, the biomass must be entered again. 
    ReEntry = 1
    while ReEntry:
        ReEntry = 0
        try:
            capacity = int(input('biomass [g/L] for the experiment: '))
            if capacity > Mutant._Mutant__BiomassMax or not capacity:
                raise ValueError('wrong input')
        except ValueError:
            print('Error, no biomass was set or unexpected value or the maximum possible biomass was exceeded. Enter a value for the biomass again.')
            ReEntry = 1
        
    
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
