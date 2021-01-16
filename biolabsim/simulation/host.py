
from typing import Any, Literal, List, Optional
from pathlib import Path
from copy import copy

from .config import METABOLIC_MODEL_DIR
from .random import pick_choice, pick_integer
from .strain import Strain, WildtypeStrain, MutatedStrain
from .common import Sequence



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

    # Maximum biomass concentration,
    max_biomass : int

    # Optimal growth temperature.
    opt_growth_temp : int

    # Factor which influences the range of the promoter strength.
    promoter_str_factor : int

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
        self.promoter_str_factor = pick_integer(30,50) # explanation see Plot_ExpressionRate
        self.opt_primer_len = pick_integer(16,28) # unit: nt, source: https://link.springer.com/article/10.1007/s10529-013-1249-8



    def clone_with_mutation ( self, name:str ) -> 'Host' :
        """
        Will return a cloned host that contains one mutation.
        """
        if self.strain is None :
            raise HostHasNoStrain()

        # Copy this host and change the mutated pieces.
        new_host = copy(self) # Shallow copy because strain is re-created.
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



    def print_biotech_setting ( self)  -> None:
        '''Report of all properties defined in the biotech experiment.'''
        print("{}: {}".format( "Host", self.name ))
        print("{}: {}".format( "Resources", self.resources ))
        print("{}: {}".format( "Substrate", self.substrate ))
        print("{}: {}".format( "Strain", [ self.strain.name if self.strain != None else "<no-strain>" ] ))



    # TODO: Code below is not integrated.



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


