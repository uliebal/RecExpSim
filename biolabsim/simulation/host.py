
from typing import Any, Literal, List, Optional
from pathlib import Path
from copy import deepcopy

from ..config import METABOLIC_MODEL_DIR
from ..random import pick_choice, pick_integer
from ..common import Sequence
from .strain import Strain, WildtypeStrain, MutatedStrain
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord



ProductionPhase = Literal[ 'exponential', 'stationary' ]



class HostHasNoStrain (Exception) :
    """Raised when a Host has no Strain but wants to use one."""
    def __init__ ( self ) :
        super().__init__("Host has no Strain but wants to perform an action with a Strain.")



class Host:
    """
    The 'Host' class stores information about the organism are placed.
    The specific strain of the host can be unknown.
    """

    # Production phase: either during growth phase or stationary phase.
    prod_phase : ProductionPhase

    # Resource available in the host.
    resources : int

    # TODO: Substrate seems to be always None.
    substrate : Any

    # Maximum biomass concentration.
    max_biomass : int

    # Optimal growth temperature.
    opt_growth_temp : int

    # Factor which influences the range of the promoter strength.
    infl_prom_streng : int

    # Optimal Primer length.
    opt_primer_len : int

    # The strain of this host. If unknown this will be set to None.
    strain : Optional[Strain]



    def __init__ ( self, name:str, max_biomass:int, strain:Optional[Strain] = None ) :

        self.name = name

        self.max_biomass = max_biomass
        self.strain = strain

        self.prod_phase = pick_choice(['exponential','stationary'])
        self.resources = 40
        self.substrate = None

        self.opt_growth_temp = pick_integer(25,40) # unit: degree celsius, source: https://application.wiley-vch.de/books/sample/3527335153_c01.pdf
        self.infl_prom_streng = pick_integer(30,50) # explanation see Plot_ExpressionRate
        self.opt_primer_len = pick_integer(16,28) # unit: nt, source: https://link.springer.com/article/10.1007/s10529-013-1249-8



    def clone_with_mutation ( self, name:str ) -> 'Host' :
        """
        Will return a cloned host that contains one mutation.
        """
        if self.strain is None :
            raise HostHasNoStrain()

        # Copy this host and change the mutated pieces.
        new_host = deepcopy(self) # Shallow copy because strain is re-created.
            # TODO: If there are problems of Dataframes being mutated in_place, change this to deep copy.
        new_host.strain = MutatedStrain( name=name, host_name=self.name, ref_strain=self.strain )

        return new_host



    def get_genome ( self ) -> Sequence :
        """
        Get the genome associated with this host. Easier access to the genome for end-users.
        """
        if self.strain is None :
            raise HostHasNoStrain()
        return self.strain.genome



    # REFACTOR: Some separate but related points:
    #   - An important concept on Object-Oriented Programming (using classes) is that
    #     each class should have a 'single responsibility'.
    #     For example, a Host should only be do things that are normal for hosts: hold information
    #     about its existence, clone itself, eat, etc. This premise of single-responsibility
    #     prevents the Host from having properties such as `var_Library` because it means that
    #     Hosts normally keep track of their lineage. Even knowing what an Experiment ID is might
    #     already be too much responsibility for a Host.
    #   - Whenever possible, prefer returning values than outputting them. Ideally, the entire
    #     `biolabsim` module should have no prints as side-effects from calculations. That way,
    #     it is the end-user's decision if they want to use the values for further calculations
    #     or if they want to print them on the screen/graph/file. This also applies for plotting.
    #     End-user in this case is whoever is using `biolabsim`: either someone importing the
    #     module or the code inside the notebook.
    #   - Methods like `Test(self)` mention a very specific type of Host. Given the specificity,
    #     this method should better be external to this class definition.
    #   - Having some accessor methods to the strain are possible to make the Host API easier to use.
    #     Like `host.get_genome()` just returning `host.strain.genome`.



    def show_BiotechSetting (self)  -> None:
        '''Report of all properties defined in the biotech experiment.'''
        print("{}: {}".format( "Host", self.name ))
        print("{}: {}".format( "Resources", self.resources ))
        print("{}: {}".format( "Substrate", self.substrate ))
        print("{}: {}".format( "Strain", self.strain.name if self.strain != None else "<no-strain>" ))

    def show_Library(self):
        '''Report of clones and their performance.'''
        for Clone_ID, Clone_info in self.var_Library.items():
            print("\nClone ID: {}".format(Clone_ID))
            for key in Clone_info:
                print('{}: {}'.format(key, Clone_info[key]))

    def show_TargetExpressionRate(self):
        '''Function to calculate the maximum possible expression rate and to tell the students what the minimum rate should be.'''
        from ..measurement.physiology import Express_Max
        achieveExpRate = Express_Max(self)
        print('Maximum possible expression: {}'.format(achieveExpRate))

    # REFACTOR: The best place to make decisions depending on the Host's identity is when
    # initializing the host, such as the case with the preset Ecol and Pput. That initialization
    # should contain all necessary info such that in the future no more 'if's have to be done to
    # decide what to do.
    # Suggestion for this case: Set the OptimalPromoterStrength on initialization and remove the
    # 'if's on all helper functions that check for the name. The Host class might need to have more
    # properties such as optimal promoter strength and expression predictor paths.

    def plot_ReferencePromoterStrength(self):
        '''Function to plot the promoter strength of the optimal sequence additionally as reference.'''
        import matplotlib.pyplot as plt

        factor = self.infl_prom_streng
        # Values see init function at the beginning
        if self.name == 'Ecol':
            OptimalPromoterStrength = round(0.057 * factor,2)
        elif self.name == 'Pput':
            OptimalPromoterStrength = round(0.04 * factor,2)
        # plot of maximum Promoter strength together with GC content
        # GC-content is the same for of both optimal sequences.
        plt.plot(0.575, OptimalPromoterStrength, marker = '*', color = 'green', markersize = 10)

    def make_TempGrowthExp(self, CultTemps, ExpID=1):
        '''Experiment to determine optimal growth rate. The experiment runs until the maximum biomass is reached.'''
        from ..measurement.physiology import Help_TempGrowthExp
        Help_TempGrowthExp(self, CultTemps, ExpID=1)

    def make_Cloning(self, Clone_ID, Promoter, Primer, Tm):
        from ..manipulation.genetic import Help_Cloning
        Help_Cloning(self, Clone_ID, Promoter, Primer, Tm)

    def measure_PromoterStrength(self, Clone_ID):
        from ..measurement.physiology import Help_MeasurePromoterStrength
        Help_MeasurePromoterStrength(self, Clone_ID)

    def measure_ProductionExperiment(self, Clone_ID, CultTemp, GrowthRate, Biomass):
        from ..measurement.physiology import Help_ProductionExperiment
        Help_ProductionExperiment(self, Clone_ID, CultTemp, GrowthRate, Biomass)

    def Test(self):
        from ..auxfun import Help_Test
        self.name = 'Ecol'
        self.resources = 40
        self.substrate = None
        self.var_Library = {}
        self.strain = None
        self.infl_prom_streng = 30 # explanation see Plot_ExpressionRate
        self.opt_growth_temp = 30 # unit: degree celsius, source: https://application.wiley-vch.de/books/sample/3527335153_c01.pdf
        self.opt_primer_len = 20 # unit: nt, source: https://link.springer.com/article/10.1007/s10529-013-1249-8
        self.max_biomass = 50
        Help_Test(self)


    def Export_Genome(self):
        '''
        Exports the whole genome sequence as fasta file.
        '''
        Genome_File = '{}_{}_Genome-Sequence.fasta'.format(self.name, self.strain.name)
        with open(Genome_File, 'w') as file:
            file.write('>FullSeq\n{}'.format(self.get_genome()))
        print('Exported genome sequence data as: {}'.format(Genome_File))

    def Export_ORFAnnotation(self):
        '''
        Saves the ORFs in fasta file.
        '''
        # Using biopython as here:
        # http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec33
        # http://python.omics.wiki/biopython/examples/read-fasta
        ORF_File = '{}_{}_ORF-Sequence.fasta'.format(self.name, self.strain.name)
        with open(ORF_File, 'w') as file:
            for Idx, Row in self.strain.genes_df.iterrows():
                mySeq = SeqRecord(Seq(Row['ORF']), id=Row['RctID'])
                # write new fasta file
                r=SeqIO.write(mySeq, file, 'fasta')
                if r!=1: print('Error while writing sequence:  ' + mySeq.id)
        print('Exported ORF sequence data as: {}'.format(ORF_File))

        