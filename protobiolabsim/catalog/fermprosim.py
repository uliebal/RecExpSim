"""
All necessary user-interactible classes that are used for the Fermentation Process Simulation.
"""

from __future__ import annotations
from typing import Optional

from ..random import set_seed, pick_uniform, pick_integer
from ..experiment import Experiment
from ..organism import Organism

from ..config import DATADIR
from ..extensions.modules.growth_behaviour import GrowthBehaviour
from ..extensions.modules.fermentation_model import (
    MonodModel, ParamSet, ConditionSet, randomize_params, randomize_cond, FermentationResult,
    OperationMode
)



class BaseMonodOrganism (Organism) :

    growth: GrowthBehaviour

    monod: MonodModel



    def __init__ (
            self, exp: Experiment, seed:Optional[int],
            opt_growth_temp: float,
            opt_growth_ph: float,
            monod_params: ParamSet,
            monod_conds: ConditionSet,
            operation_mode: OperationMode ) :
        super().__init__( exp=exp, seed=seed )

        self.growth = GrowthBehaviour(
            org=self,
            opt_growth_temp=opt_growth_temp,
            opt_growth_ph=opt_growth_ph,
            max_biomass=1,  # NOTE: Max_Biomass is not used by MonodModel.
        )

        self.monod = MonodModel(
            org=self,
            growth=self.growth,
            params=monod_params,
            conditions=monod_conds,
            operation_mode=operation_mode,
        )



    def clone ( self ) -> BaseMonodOrganism :
        """ Clone the organism while using the same randomized values. """
        return BaseMonodOrganism(
            exp=self.exp,
            opt_growth_temp = self.growth.opt_growth_temp,
            opt_growth_ph = self.growth.opt_growth_ph,
            monod_params = self.monod.params,
            monod_conds = self.monod.conditions
        )



    def print_status ( self ) -> None :
        print("Organism Information:")
        print("  opt_growth_temp = {}".format( self.growth.opt_growth_temp ))
        print("  opt_growth_ph = {}".format( self.growth.opt_growth_ph ))
        print("  monod params = {}".format( self.monod.params ))
        print("  monod conds = {}".format( self.monod.conditions ))



    def calc_monod_kinetics ( self ) -> FermentationResult :
        result = self.monod.calculate_monod()
        # NOTE: This result is given to the user. If we want to cache the result,
        # then this would be the place to do it.
        return result



class SomeAcidicOrganism (BaseMonodOrganism) :
    def __init__ ( self, exp: Experiment, seed: Optional[int] = None ) :
        set_seed(seed)
        super().__init__(
            exp=exp, seed=seed,
            opt_growth_temp=30, # NOTE: Example: pick_integer(25,40),
            opt_growth_ph=5, # NOTE: Could be randomized like temp.
            monod_params=randomize_params(seed=seed),
            monod_conds=randomize_cond(seed=seed, duration=20),
            operation_mode='batch'
        )



class SomeColdOrganism (BaseMonodOrganism) :
    def __init__ ( self, exp: Experiment, seed: Optional[int] = None ) :
        set_seed(seed)
        super().__init__(
            exp=exp, seed=seed,
            opt_growth_temp=21, # NOTE: Example: pick_integer(25,40),
            opt_growth_ph=6, # NOTE: Could be randomized like temp.
            monod_params=randomize_params(seed=seed),
            monod_conds=randomize_cond(seed=seed, duration=20),
            operation_mode='batch'
        )



class FermExperiment (Experiment) :

    def __init__ ( self ) :
        super().__init__()

    def create_acidic_organism ( self, seed:Optional[int] ) :
        return SomeAcidicOrganism( exp=self, seed=seed )

    def create_cold_organism ( self, seed:Optional[int] ) :
        return SomeColdOrganism( exp=self, seed=seed )
