
from typing import Any, Literal, List
from pathlib import Path

from .config import MODEL_DIR
from .random import pick_choice, pick_integer
from .strain import Strain, WildtypeStrain, MutatedStrain


ProductionPhase = Literal[ 'exponential', 'stationary' ]




class Host:
    """
    The 'Host' class stores all information about the organism and the integrated recombinant protein.

    Multiple mutations can be stacked to form a mutated Host. The Host will keep track of all
    Strains.
    """

    # Production phase: either during growth phase or stationary phase.
    # TODO: This is currently never used.
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

    # A stack of strains, where the bottom strain is the wildtype and all other strains are
    # mutations on top of the previous strains.
    # TODO: Does it make sense for a host to have multiple Strains ever. Or is the Host the Strain?
    strains : List[Strain]



    def __init__ ( self, name:str, max_biomass:int, metabolic_model_path:Path ) :

        self.name = name

        self.prod_phase = pick_choice(['exponential','stationary'])
        self.resources = 40
        self.substrate = None

        # TODO: These variables are mostly used in BioLabSimFun
        self.opt_growth_temp = pick_integer(25,40) # unit: degree celsius, source: https://application.wiley-vch.de/books/sample/3527335153_c01.pdf
        self.promoter_str_factor = pick_integer(30,50) # explanation see Plot_ExpressionRate
        self.opt_primer_len = pick_integer(16,28) # unit: nt, source: https://link.springer.com/article/10.1007/s10529-013-1249-8

        # Initiating metabolic network by adding a default WT genome-metabolism connection.
        self.strains = [
            WildtypeStrain( name="WT", host=self, model_path=metabolic_model_path )
        ]



    def wt_strain ( self ) -> Strain :
        """
        Returns the first strain that is always the WT strain.
        """
        return self.strains[0]



    def last_strain ( self ) -> Strain :
        """
        Returns the top-most strain that holds all applied mutations on the wildtype
        Most of the times the last mutation, but sometimes also the wildtype itself.
        """
        return self.strains[-1]



    def add_mutation ( self, name:str ) -> None :
        '''
        Adding mutation on top of the current strain stack.
        '''
        new_strain = MutatedStrain( name=name, host=self, ref_strain=self.last_strain() )
        self.strains.append( new_strain )



    def print_biotech_setting ( self)  -> None:
        '''Report of all properties defined in the biotech experiment.'''
        print("{}: {}".format( "Host", self.name ))
        print("{}: {}".format( "Resources", self.resources ))
        print("{}: {}".format( "Substrate", self.substrate ))
        print("{}: {}".format( "StrainLibrary", [ s.name for s in self.strains ] ))



    # TODO: Not tested or sure if is correct.
    def show_Library (self) -> None :
        '''Report of clones and their performance.'''
        for clone in self.strains.items():
            Clone_ID = clone.name
            print("\nClone ID: {}".format(Clone_ID))
            for key in Clone_info:
                print('{}: {}'.format(key, Clone_info[key]))

