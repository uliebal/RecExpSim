class Host:
    '''The 'Host' class stores all information about the organism and the integrated recombinant protein.''' 
    
#     from BioLabSimFun import SequenceRandomizer_Single
    from random import randint
    # random assignment of the production phase, either during growth phase or stationary phase
    __ProdPhase = 'exponential' if randint(0,1)==0 else 'stationary'
    # resources, e.g. money, for conducting tests
    __Resources = 40
    __BiomassMax = None
    
    def __init__(self, Host, Metabolism=False):
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
        # adding default genome-metabolism connection called WT
        if Metabolism:
            Name = 'WT'
            myWT = Strain(Name)
            myWT.make_WT(self)
            setattr(self, Name, myWT)
            self.var_StrainLibrary = [Name]
        
    def add_Strain(self, Name):
        '''
        Adding strains with different genome to the host.
        '''
        myStrain = Strain(Name)
        myStrain.make_Mutant(self, NumberEnzymes=3, Target='-10', NumberMutations=2)
        setattr(self, Name, myStrain)
        self.var_StrainLibrary.append(Name)

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
                
    def show_TargetExpressionRate(self):
        '''Function to calculate the maximum possible expression rate and to tell the students what the minimum rate should be.'''
        from BioLabSim.ModuleMeasureOrganism.Physiology import Express_Max
        achieveExpRate = Express_Max(self)
        print('Maximum possible expression: {}'.format(achieveExpRate))
            
    def make_TempGrowthExp(self, CultTemps, ExpID=1):
        '''Experiment to determine optimal growth rate. The experiment runs until the maximum biomass is reached.'''
        from BioLabSim.ModuleMeasureOrganism.Physiology import Help_TempGrowthExp
        Help_TempGrowthExp(self, CultTemps, ExpID=1)

    def make_Cloning(self, Clone_ID, Promoter, Primer, Tm):
        from BioLabSim.ModuleManipulateOrganism.GeneticChanges import Help_Cloning        
        Help_Cloning(self, Clone_ID, Promoter, Primer, Tm)
        
    def measure_PromoterStrength(self, Clone_ID):
        from BioLabSim.ModuleMeasureOrganism.Physiology import Help_MeasurePromoterStrength
        Help_MeasurePromoterStrength(self, Clone_ID)
        
    def measure_ProductionExperiment(self, Clone_ID, CultTemp, GrowthRate, Biomass):
        from BioLabSim.ModuleMeasureOrganism.Physiology import Help_ProductionExperiment
        Help_ProductionExperiment(self, Clone_ID, CultTemp, GrowthRate, Biomass)
        
    def Test(self):
        from BioLabSim.AuxFun import Help_Test
        self.var_Host = 'Ecol'
        self.var_Resources = 40
        self.var_Substrate = None
        self.var_Library = {}
        self.info_Strain = {}
        self.__InflProStreng = 30 # explanation see Plot_ExpressionRate 
        self.__OptTemp = 30 # unit: degree celsius, source: https://application.wiley-vch.de/books/sample/3527335153_c01.pdf
        self.__OptPrLen = 20 # unit: nt, source: https://link.springer.com/article/10.1007/s10529-013-1249-8
        self.__BiomassMax = 50
        Help_Test(self)



class Strain:
    '''
    The 'strain' class stores information of genome and other associated metabolic information
    '''
    def __init__(self, Name):
        import os
        from BioLabSim.ModuleVirtualOrganism.Metabolism import Help_LoadCobra
        
        print('Initiating metabolic network')

        ModelPath = os.path.join('Models','e_coli_core.xml')

        self.ID = Name

        self.var_Model = Help_LoadCobra(Path = ModelPath)
        
    def make_WT(self, Host, GenomeSize = 500, GCcont = .6):
        from BioLabSim.ModuleVirtualOrganism.Metabolism import Help_GeneAnnotator, Help_Expr2Flux, Help_FluxCalculator
        from BioLabSim.ModuleVirtualOrganism.Genome import Help_GenomeGenerator
        
        self.var_Host = getattr(Host, 'var_Host')
        self.info_GenesDF = Help_GeneAnnotator(self.var_Host, self.var_Model) # self.__GenesDF = 
        self.info_GenesDF['Fluxes'], self.info_Objective = Help_FluxCalculator(self, Host)
        self.info_GenesDF['Expr2Flux'] = Help_Expr2Flux(self.info_GenesDF) # self.__GenesDF = 
        self.info_Genome = Help_GenomeGenerator(self.info_GenesDF, GenomeSize, GCcont)
        self.var_GenomeSize = len(self.info_Genome)
        self.var_GCcont = round((self.info_Genome.count('G') + self.info_Genome.count('C'))/self.var_GenomeSize,2)
        
    def make_Mutant(self, Host, NumberEnzymes=3, Target='-10', NumberMutations=2):
        from BioLabSim.ModuleVirtualOrganism.Metabolism import Help_Expr2Flux, Help_StrainCharacterizer, Help_FluxCalculator
        from BioLabSim.ModuleManipulateOrganism.GeneticChanges import Help_MutActProm
        
        RefGenome = Host.WT.info_Genome
        RefGenDF = Host.WT.info_GenesDF
        Host_strain = Host.var_Host
        GenomeMut, GenDFMut = Help_MutActProm(RefGenome, RefGenDF, NumberEnzymes=3, Target='-10', NumberMutations=2)
        self.info_GenesDF = GenDFMut
        self.info_Genome = GenomeMut
#         RctNewDF = Help_StrainCharacterizer(Host_strain, RefGenDF, RefGenome, GenomeMut, self.var_Model)
#         self.info_GenesDF['Fluxes'] = Help_FluxCalculator(self.var_Model, RctNewDF)
        self.info_GenesDF['Fluxes'], self.info_Objective = Help_FluxCalculator(self, Host, True)
        self.info_GenesDF['Expr2Flux'] = Help_Expr2Flux(self.info_GenesDF)
        self.var_GenomeSize = len(self.info_Genome)
        self.var_GCcont = round((self.info_Genome.count('G') + self.info_Genome.count('C'))/self.var_GenomeSize,2)


        
        



