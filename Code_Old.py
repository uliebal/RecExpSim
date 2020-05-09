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
