"""
GenomeLibrary is a module that stores the multiple genes an Organism may have.

However, it only stores the references to those Genes. The actual storage place of the genes is
in `registry.gene` which has special ways to allow for copy-on-write.
"""

from __future__ import annotations
from copy import copy
from typing import Optional
from random import randint
from typing import TypedDict, Tuple, List
from collections import namedtuple

import numpy as np
import pandas as pd

from ...random import pick_uniform, pick_normal
from ...organism import Organism
from ...module import Module
#from .growth_record import Growth



class GrowthBehaviour ( Module ) :

    opt_growth_temp: int
    max_biomass: int


    # infl_prom_streng = randint(30,50) # explanation see Plot_ExpressionRate
    # # optimal growth temperature, randomly assigned
    # opt_growth_temp = randint(25,40) # unit: degree celsius, source: https://application.wiley-vch.de/books/sample/3527335153_c01.pdf
    # # optimal Primer length, randomly assigned
    # opt_primer_len = randint(16,28) # unit: nt, source: https://link.springer.com/article/10.1007/s10529-013-1249-8
    # maximum biomass concentration, the limits for Ecol were set as shown below and the values for Pput were adjusted according to the ratio of the maximum promoter strengths (0.057/0.04) of the optimal sequences (see expression measurement issue).
    # if name == 'Ecol':
    #     max_biomass = randint(30,100) # unit: in gDCW/l, source (german): https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=2&cad=rja&uact=8&ved=2ahUKEwjzt_aJ9pzpAhWGiqQKHb1jC6MQFjABegQIAhAB&url=https%3A%2F%2Fwww.repo.uni-hannover.de%2Fbitstream%2Fhandle%2F123456789%2F3512%2FDissertation.pdf%3Fsequence%3D1&usg=AOvVaw2XfGH11P9gK2F2B63mY4IM
    # elif name == 'Pput':
    #     max_biomass = randint(45,145) # unit: in gDCW/l, source 1: https://onlinelibrary.wiley.com/doi/pdf/10.1002/bit.25474, source 2: https://link.springer.com/article/10.1385/ABAB:119:1:51


    # TODO: I'm not happy about this constructor approach with ref. I don't like that all parameters need
    # to be optional just in case a 'ref' is passed. Maybe do clone(dep_mods,ref) which itself
    # calls the init with all parameters passed by reference. And then deconstruct the 'params' to
    # use kwargs.
    def __init__ (
        self, org:Organism, ref:Optional[GrowthBehaviour] = None,
        opt_growth_temp:Optional[int] = None, max_biomass:Optional[int] = None
    ) :
        super().__init__(org,ref)

        if ref is not None :
            self.opt_growth_temp = ref.opt_growth_temp
            self.max_biomass = ref.max_biomass
        else :
            self.opt_growth_temp = opt_growth_temp # TODO: no validation
            self.max_biomass = max_biomass # TODO: no validation



    def grow ( self, CultTemps:list[int] ) -> Tuple[ pd.DataFrame, list ]:
        """
        TODO: Ideally, growth should return a table with enough information to contain the
        biomass results and to infer the loading times. Loading times is a synthetic construct
        that should be contained in the Organism, while the Module should only contain efficient
        calculations. The Organism is an API for the User and may add interaction such as
        Resource management and loading times.
        """
        CultTemps = np.array(CultTemps)

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

        PauseEntry = namedtuple('PauseEntry', 'exp loading_len')
        pauses = []

        #computing of biomass data and updating of DataFrame
        for i in range(len(CultTemps)):

            if pick_uniform(0,1) > 0.1: # in 10% of cases the culture does not grow (failure of the experiment)
                r = Help_GrowthConstant(OptTemp, CultTemps[i])
                # the result can reach very small values, which poses downstream problems, hence the lowest value is set to 0.05
                if r > 0.05: # under process conditions it might be realistic, source : https://www.thieme-connect.de/products/ebooks/pdf/10.1055/b-0034-10021.pdf
                    duration = d_mult * 1/r * np.log((capacity - P0)/P0) + 1
                else:
                    duration = 7
                t = np.arange(np.minimum(73, duration))

                # biomass data is calculated according to https://en.wikipedia.org/wiki/Logistic_function
                mu = capacity / (1 + (capacity-P0) / P0 * np.exp(-r * t))
                sigma = 0.1*mu

                exp_TempGrowthExp = [pick_normal(mu[k], sigma[k]) for k in range(len(mu))]

                loading_len = len(t)
                exp = ' of exp.{} at {} °C'.format(i+1, (CultTemps[i]))
                pauses.append(PauseEntry( exp, loading_len ))

            else:
                mu = P0
                sigma = 0.08*mu
                exp_TempGrowthExp = [pick_normal(mu, sigma) for i in range(7)] # if cells haven't grown, the measurement is only continued for 6h

                loading_len = 7
                exp = ' of exp.{} at {} °C'.format(i+1, (CultTemps[i]))
                pauses.append(PauseEntry( exp, loading_len ))


            new_df = pd.DataFrame({'exp.{} biomass conc. at {} °C'.format(i+1, (CultTemps[i])): exp_TempGrowthExp})
            df.update(new_df)

        return ( df, pauses )










def Help_GrowthConstant(OptTemp, CultTemp, var=5):
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
#     if Biomass > Mutant._Mutantmax_biomass or not Biomass:
#         print('Error, no biomass was set or unexpected value or the maximum possible biomass was exceeded. Enter a value for the biomass again.')

    # Equation for calculating the maximum slope
    # https://www.tjmahr.com/anatomy-of-a-logistic-growth-curve/
    GrowthMax = Biomass * growth_rate_const / 4

    return GrowthMax