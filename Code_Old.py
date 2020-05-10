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

    #def show_Primer_DeviationOptimalLength(self, PrLen):
        #import numpy as np
        #OptLen = self._Mutant__OptPrLen
        #Devi = (np.absolute(OptLen - PrLen)/OptLen)*100
        #print('The deviation from the optimum length is {} %.'.format(Devi.round(2)))
        #if PrLen <= 30:
            #print('The length of the primer does not exceed 30 nt.')
        #else:
            #print('The primer is too long (>30 nt).')

#         [Clone_ID]['Promoter_GC-content'] = (Promoter.count('C') + Promoter.count('G')) / len(Promoter)
#  = {Clone_ID: {'Promoter_GC-content': (Promoter.count('C') + Promoter.count('G')) / len(Promoter)}}
#         self.var_Promoter = Promoter
#         self.var_GC_content = (self.var_Promoter.count('C') + self.var_Promoter.count('G')) / len(self.var_Promoter)

#     def __add_RandomPromoter(self):
#         self.var_Promoter = SequenceRandomizer_Single()
#         self.var_GCcontent = (self.var_Promoter.count('C') + self.var_Promoter.count('G')) / len(self.var_Promoter)

#         if draw_plot:
#             WaitingTime = 30/r
#             for i in range (len(t)):
#                 pl.clf()
#                 pl.figure(figsize = (5,3), dpi = 120)
#                 pl.xlim(0, t[i])
#                 pl.ylim(0, 1.1*np.max(exp_TempGrowthExp))
#                 pl.xlabel("time [h]")
#                 pl.ylabel("biomass concentration [g/L]")
#                 pl.plot(t[0:i], exp_TempGrowthExp[0:i], linestyle = '--') 
#                 display.display(pl.gcf())
#                 display.clear_output(wait=True)
#                 #pl.pause(1)
#                 time.sleep(WaitingTime/n) # total waiting time independent of n  
        
        #Help_ExportToExcel(self, 'Strain_characterization', 'different Temp',
                  #t, exp_TempGrowthExp, 'time', f'biomass {CultTemp}') 
