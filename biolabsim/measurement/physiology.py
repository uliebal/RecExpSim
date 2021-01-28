def Help_TempGrowthExp(Host, CultTemps, ExpID=1):
    '''Experiment to determine optimal growth rate. The experiment runs until the maximum biomass is reached.'''
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import time

    from ..random import pick_uniform, pick_normal
    from ..auxfun import Help_Progressbar, Error_Resources

    if Host.resources > 0:

        capacity = Host.max_biomass
        # the time of the half maximum population (inflection point) is calculated according to here:
        # https://opentextbc.ca/calculusv2openstax/chapter/the-logistic-equation/
        d_mult = 2 # we multiply the inflection point with 'd_mult' to increase cultivation time
        P0 = 0.1

        # determine time vector with maximum length:
        OptTemp = Host.opt_growth_temp
        # Selecting the temperature that is most distant from the optimal temperature
        Temp_tmax = CultTemps[np.argmax(np.absolute(CultTemps-OptTemp))]
        # using the worst temperature to calculate lowest growth rate
        r_tmax = Help_GrowthConstant(Host, Temp_tmax)
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
            if Host.resources > 0:
                wait = 0.01 # has to be adjusted, waiting time for loading bar

                if pick_uniform(0,1) > 0.1: # in 10% of cases the culture does not grow (failure of the experiment)
                    r = Help_GrowthConstant(Host, CultTemps[i])
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

                    loading_time = wait * len(t)
                    exp = ' of exp.{} at {} °C'.format(i+1, (CultTemps[i]))
                    Help_Progressbar(45, loading_time, exp)

                else:
                    mu = P0
                    sigma = 0.08*mu
                    exp_TempGrowthExp = [pick_normal(mu, sigma) for i in range(7)] # if cells haven't grown, the measurement is only continued for 6h

                    loading_time = wait * 7
                    exp = ' of exp.{} at {} °C'.format(i+1, (CultTemps[i]))
                    Help_Progressbar(45, loading_time, exp)

            else:
                Error_Resources()
                return

            new_df = pd.DataFrame({'exp.{} biomass conc. at {} °C'.format(i+1, (CultTemps[i])): exp_TempGrowthExp})
            df.update(new_df)

            Host.resources -= 1

        # Export DataFrame to Excel
        with pd.ExcelWriter('Tst_Strain_characterization_{}.xlsx'.format(ExpID)) as writer:
            df.to_excel(writer, sheet_name='different temp')

    print('Tst_Strain_characterization_{}.xlsx'.format(ExpID))
    else:
        Error_Resources()
        return


def Help_MeasurePromoterStrength(Host, Clone_ID):
    from .sequencing import Help_PromoterStrength
    if Host.resources > 0:
        if hasattr(Host, 'var_Library'):
            if Clone_ID in Host.var_Library:
                factor = Host.infl_prom_streng
                Sequence = Host.var_Library[Clone_ID]['Promoter_Sequence']
                Host.var_Library[Clone_ID]['Promoter_Strength'] = round(Help_PromoterStrength(Host.var_Host, Sequence) * factor, 2)
                Host.resources -= 1
            else:
                print('Error, Clone ID does not exist. Choose existing Clone ID.')
        else:
            print('Error, no promoter library available. Perform a cloning first.')
    else:
        Error_Resources()

def Help_ProductionExperiment(Host, Clone_ID, CultTemp, GrowthRate, Biomass, accuracy_Test=.9):
    import numpy as np
    if Host.resources > 2: # three resources will be deducted
        # the final experiment can only be performed after at least one sequence has been cloned and tested:
        if hasattr(Host, 'var_Library'):
            if Clone_ID in Host.var_Library:
                # testing whether the determined maximum biomass and the determined maximum growth rate are close to the actual ones
                Tst_Biomass = 1 - np.abs(Biomass-Host.max_biomass) / Host.max_biomass
                Tst_Temp = 1-np.abs(GrowthRate-Help_GrowthConstant(Host, Host.opt_growth_temp))/Help_GrowthConstant(Host, Host.opt_growth_temp)
                if  Tst_Biomass > accuracy_Test and Tst_Temp > accuracy_Test:
                    # Growth rate was only checked, for the calculation the rate resulting from the temperature is used
                    r = Help_GrowthConstant(Host, CultTemp)
                    GrowthMax = Growth_Maxrate(r, Biomass)
                    Host.var_Library[Clone_ID]['Expression_Temperature'] = CultTemp
                    Host.var_Library[Clone_ID]['Expression_Biomass'] = Biomass
                    Host.var_Library[Clone_ID]['Expression_Rate'] = round(GrowthMax * Host.var_Library[Clone_ID]['Promoter_Strength'],2)
                    Host.resources -= 3
                else:
                    print('Maximum biomass and/or maximum growth rate are incorrect.')
            else:
                print('Error, Clone ID does not exist. Choose existing Clone ID.')
        else:
            print('Error, no promoter sequence has been cloned and tested yet. Perform a cloning first and then test the expression with "Make_MeasurePromoterStrength(Clone_ID)".')

    else:
        Error_Resources()

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

    OptTemp = Host.opt_growth_temp
    r_pdf = norm(OptTemp, var)
    # calculation of the growth rate constant, by picking the activity from a normal distribution

    growth_rate_const = r_pdf.pdf(CultTemp) / r_pdf.pdf(OptTemp)

    return growth_rate_const



def Growth_Maxrate(growth_rate_const, Biomass):
    '''The function calculates the maximum slope during growth.
        Arguments:
            growth_rate_const: float, maximum growth rate constant
            Biomass:           float, maximum biomass concentration as carrying capacity
        Output:
            growth_rate_max:   float, maximum growth rate
    '''
    # biomass checks
#     if Biomass > Mutant._Mutantmax_biomass or not Biomass:
#         print('Error, no biomass was set or unexpected value or the maximum possible biomass was exceeded. Enter a value for the biomass again.')

    # Equation for calculating the maximum slope
    # https://www.tjmahr.com/anatomy-of-a-logistic-growth-curve/
    GrowthMax = Biomass * growth_rate_const / 4

    return GrowthMax

def Express_Max(Host):
    '''
    Determines maximum possible expression values
    '''
    BiomassMax = Host.max_biomass
    OptTemp = Host.opt_growth_temp
    factor = Host.infl_prom_streng
    # Values see init function at the beginning
    if Host.var_Host == 'Ecol':
        MaximumPromoterStrength = round(0.057 * factor,2)
    elif Host.var_Host == 'Pput':
        MaximumPromoterStrength = round(0.04 * factor,2)
    r = Help_GrowthConstant(Host, OptTemp)
    GrowthMax = Growth_Maxrate(r, BiomassMax)
    achieveExpRate = round(GrowthMax * MaximumPromoterStrength,2)

    return achieveExpRate