"""
TODO: GrowthBehaviour needs to be reviewed.
"""

from __future__ import annotations
from copy import copy
from protobiolabsim.experiment import Experiment
from typing import Optional
from typing import TypedDict, Tuple, List
from collections import namedtuple
from dataclasses import dataclass

import numpy as np
import pandas as pd

from ...common import Outcome
from ...random import pick_uniform, pick_normal
from ...organism import Organism
from ...module import Module
from ..records.gene.gene import Gene
from ..utils import Help_GrowthConstant, Growth_Maxrate
from .genome_expression import GenomeExpression
#from .growth_record import Growth



class GrowthBehaviour ( Module ) :

    # Dependent module Genome Expression holding the genes to express.
    genexpr: GenomeExpression

    opt_growth_temp: int
    max_biomass: int


    # TODO: I'm not happy about this constructor approach with ref. I don't like that all parameters need
    # to be optional just in case a 'ref' is passed. Maybe do clone(dep_mods,ref) which itself
    # calls the init with all parameters passed by reference. And then deconstruct the 'params' to
    # use kwargs.
    def __init__ (
        self, org:Organism, genexpr:GenomeExpression,
        opt_growth_temp:int, max_biomass:int
    ) :
        super().__init__(org)

        self.genexpr = genexpr

        self.opt_growth_temp = opt_growth_temp # TODO: no validation
        self.max_biomass = max_biomass


    def clone ( self, org:Organism ) -> GrowthBehaviour :
        return GrowthBehaviour(
            org=org,
            opt_growth_temp= self.opt_growth_temp,
            max_biomass= self.max_biomass,
        )



    def Make_TempGrowthExp ( self, CultTemps:list[int] ) -> Tuple[ pd.DataFrame, list ]:
        """
        TODO: Ideally, growth should return a table with enough information to contain the
        biomass results and to infer the loading times. Loading times is a synthetic construct
        that should be contained in the Organism, while the Module should only contain efficient
        calculations. The Organism is an API for the User and may add interaction such as
        Resource management and loading times.
        """
        CultTemps = np.array(CultTemps)
        Exp_Duration = 48

        capacity = self.max_biomass
        # the time of the half maximum population (inflection point) is calculated according to here:
        # https://opentextbc.ca/calculusv2openstax/chapter/the-logistic-equation/
        d_mult = 2 # we multiply the inflection point with 'd_mult' to increase cultivation time
        P0 = 0.1

        # determine time vector with maximum length:
        OptTemp = self.opt_growth_temp
        # Selecting the temperature that is most distant from the optimal temperature
        Temp_tmax = CultTemps[np.argmax(np.absolute(CultTemps-OptTemp))]
        # using the worst temperature to calculate lowest growth rate
        r_tmax = Help_GrowthConstant(OptTemp, Temp_tmax)
        # using the worst temperature growth rate to compute longest simulation time, maximum set to 72 h
        duration_tmax = d_mult * 1/r_tmax * np.log((capacity - P0)/P0) + 1
        t_max = np.arange(np.minimum(Exp_Duration, duration_tmax))

        # create an empty DataFrame with t_max as first column
        col = []
        col.append('time [h]')
        for i in range(len(CultTemps)):
            col.append('exp.{} biomass conc. at {} °C'.format(i, (CultTemps[i])))

        df = pd.DataFrame(np.empty(shape=(len(t_max), len(CultTemps)+1), dtype=float), columns = col)
        df[:len(t_max)] = np.nan
        new_df = pd.DataFrame({'time [h]': t_max})
        df.update(new_df)

        PauseEntry = namedtuple('PauseEntry', 'exp loading_len')
        pauses = []

        #computing of biomass data and updating of DataFrame
        for i in range(len(CultTemps)):

            if pick_uniform(0,1) > self.org.exp.suc_rate: # TODO: Too nested. # experiment failure depending on investment to equipment
                r = Help_GrowthConstant(OptTemp, CultTemps[i])
                # the result can reach very small values, which poses downstream problems, hence the lowest value is set to 0.05
                if r > 0.05: # under process conditions it might be realistic, source : https://www.thieme-connect.de/products/ebooks/pdf/10.1055/b-0034-10021.pdf
                    duration = Exp_Duration #d_mult * 1/r * np.log((capacity - P0)/P0) + 1
                else:
                    duration = Exp_Duration
                t = np.arange(np.minimum(Exp_Duration, duration))

                # biomass data is calculated according to https://en.wikipedia.org/wiki/Logistic_function
                mu = capacity / (1 + (capacity-P0) / P0 * np.exp(-r * t))
                sigma = 0.1*mu

                exp_TempGrowthExp = [pick_normal(mu[k], sigma[k]) for k in range(len(mu))]

                loading_len = len(t)
                exp = ' of exp.{} at {} °C'.format(i, (CultTemps[i]))
                pauses.append(PauseEntry( exp, loading_len ))

            else:
                mu = P0
                sigma = 0.08*mu
                exp_TempGrowthExp = [pick_normal(mu, sigma) for i in range(Exp_Duration)] # if cells haven't grown, the measurement is only continued for 6h

                loading_len = 7
                exp = ' of exp.{} at {} °C'.format(i, (CultTemps[i]))
                pauses.append(PauseEntry( exp, loading_len ))


            new_df = pd.DataFrame({'exp.{} biomass conc. at {} °C'.format(i, (CultTemps[i])): exp_TempGrowthExp})
            df.update(new_df)

        return ( df, pauses )



    def Make_ProductionExperiment ( self, gene:Gene, CultTemp, GrowthRate, Biomass, ref_prom:str, accuracy_Test=.9 ) -> Outcome :
        """
        Return
        ======
        [ expression_rate, error_text ]
        """
        growth_const = Help_GrowthConstant(self.opt_growth_temp, self.opt_growth_temp)

        # testing whether the determined maximum biomass and the determined maximum growth rate are close to the actual ones
        if not (
            1 - np.abs(Biomass-self.max_biomass) / self.max_biomass > accuracy_Test
            and 1 - np.abs(GrowthRate-growth_const) / growth_const > accuracy_Test
        ) :
            return Outcome( 0, 'Maximum biomass and/or maximum growth rate are incorrect.' )

        # Growth rate was only checked, for the calculation the rate resulting from the temperature is used
        r = Help_GrowthConstant(self, CultTemp)
        GrowthMax = Growth_Maxrate(r, Biomass)
        PromStrength = self.genexpr.calc_prom_str( gene=gene, ref_prom=ref_from )
        AbsRate = round(GrowthMax * PromStrength,2)
        FinalRelRate = round(AbsRate/Calc_MaxExpress(self),2)
        return Outcome( FinalRelRate )



    def Calc_MaxExpress (self) -> float :
        '''Function to calculate the maximum possible expression rate.'''
        BiomassMax = self.max_biomass
        OptTemp = self.opt_growth_temp
        factor = self.genexpr.infl_prom_streng
        species_str = self.genexpr.species_prom_streng
        MaximumPromoterStrength = round(0.057 * factor,2)
        r = Help_GrowthConstant(self, OptTemp)
        GrowthMax = Growth_Maxrate(r, BiomassMax)
        return round(GrowthMax * MaximumPromoterStrength,2)