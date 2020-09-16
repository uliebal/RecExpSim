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
            self.WT = Strain()
            self.var_StrainLibrary = ['WT']
        
    def add_Strain(self,Name):
        '''
        Adding strains with different genome to the host.
        '''
        self.Name = Strain()

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





class Strain:
    '''
    The 'strain' class stores information of genome and other associated metabolic information
    '''
    def __init__(self, Host = 'Ecol', GenomeSize = 500, GCcont = .6):
        import os
        from BioLabSim.ModuleVirtualOrganism.Metabolism import Help_LoadCobra, Help_GeneAnnotator, Help_Expr2Flux
        from BioLabSim.ModuleVirtualOrganism.Genome import Help_GenomeGenerator
        print('Initiating metabolic network')
        ModelPath = os.path.join('Models','e_coli_core.xml')
        self.var_Model = Help_LoadCobra(Path = ModelPath)
        self.info_GenesDF = Help_GeneAnnotator(Host, self.var_Model) # self.__GenesDF = 
        self.info_GenesDF['Expr2Flux'] = Help_Expr2Flux(self.info_GenesDF) # self.__GenesDF = 
        self.info_Genome = Help_GenomeGenerator(self.info_GenesDF, GenomeSize, GCcont)
        self.var_GenomeSize = len(self.info_Genome)
        self.var_GCcont = round((self.info_Genome.count('G') + self.info_Genome.count('C'))/self.var_GenomeSize,2)
        

        
        
        




