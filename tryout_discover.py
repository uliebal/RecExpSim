
from biolabsim.host import Host, Strain
from biolabsim.metabolism import Help_StrainCharacterizer, measure_EnzymeLevel1
from biolabsim.genome import measure_BaseCompare
from biolabsim.organism.ecol import Ecol


myHost = Ecol()
myHost.print_biotech_setting()
# Initiating metabolic network
# Host: Pput
# Resources: 40
# Substrate: None
# Library: {}
# StrainLibrary: ['WT']


myHost.add_mutation('Mut1')
myHost.print_biotech_setting()
# Initiating metabolic network
# resetting boundaries
# Host: Pput
# Resources: 40
# Substrate: None
# Library: {}
# StrainLibrary: ['WT', 'Mut1']


myHostWT = myHost.wt_strain()
myHostWT.objective
# 0.8739215069684303


RctNewDF = Help_StrainCharacterizer(myHost.name, myHostWT.genes_df, myHostWT.genome, myHost.Mut1.genome, myHostWT.model)
RctNew = RctNewDF[RctNewDF['RctFlag']==True].index.values
print(RctNew)
# For reactions with reduced expression and positive flux: reduce the upper limit,
# For reactions with increased expression and positive flux: increase the lower limit
# for reactions with negative flux the limits are exchanged.
FluxPos = RctNew[tuple([RctNewDF.loc[RctNew, 'RefFlux']>0])]
FluxNeg = RctNew[tuple([RctNewDF.loc[RctNew, 'RefFlux']<0])]
# Finding increased and decreased fluxes
FluxInc = RctNew[RctNewDF.loc[RctNew, 'NewExpr'].values / RctNewDF.loc[RctNew, 'RefExpr'].values>1]
FluxDec = RctNew[RctNewDF.loc[RctNew, 'NewExpr'].values / RctNewDF.loc[RctNew, 'RefExpr'].values<1]
print('Positive: {}, Negative: {}'.format(FluxPos,FluxNeg))
print('Increase: {}, Decrease: {}'.format(FluxInc,FluxDec))
# [ 7 13 93]
# Positive: [13 93], Negative: [7]
# Increase: [93], Decrease: [ 7 13]


myDiff = measure_BaseCompare(myHostWT.genome, myHost.Mut1.genome)
a,b,c = measure_EnzymeLevel1(myHost, myHost.Mut1)
print(myDiff)
print(a,b,c)
# [[461, ('T', 'A')], [464, ('T', 'C')], [814, ('T', 'A')], [817, ('T', 'A')], [5636, ('A', 'T')], [5638, ('C', 'G')]]
#     RctFlag    RctID  RefExpr  NewExpr    RefFlux  Expr2Flux
# 0     False      PFK    1.059    1.059   7.477382   7.060795
# 1     False      PFL    1.726    1.726   0.000000   0.000000
# 2     False      PGI    0.607    0.607   4.860861   8.008008
# 3     False      PGK    0.482    0.482 -16.023526 -33.243830
# 4     False      PGL    1.458    1.458   4.959985   3.401910
# ..      ...      ...      ...      ...        ...        ...
# 90    False   NADH16    1.509    1.509  38.534610  25.536521
# 91    False  NADTRHD    2.523    2.523   0.000000   0.000000
# 92    False     NH4t    1.035    1.035   4.765319   4.604173
# 93     True      O2t    0.555    0.564  21.799493  39.278365
# 94    False      PDH    0.500    0.500   9.282533  18.565065

# [95 rows x 6 columns] {'lower': array([93,  7]), 'upper': array([13])} {'increase': array([93]), 'decrease': array([13,  7])}


myHostWT.model
# Name	                e_coli_core
# Memory address	    0x07fe716ead460
# Number of metabolites	72
# Number of reactions	95
# Number of groups	    0
# Objective expression	1.0*BIOMASS_Ecoli_core_w_GAM - 1.0*BIOMASS_Ecoli_core_w_GAM_reverse_712e5
# Compartments	        extracellular space, cytosol


FluxDiff = myHost.Mut1.genes_df['Fluxes']/myHostWT.genes_df['Fluxes']
np.nan_to_num(FluxDiff, nan=1, posinf=1, neginf=1)
# array([0.93983279, 1.        , 0.71771364, 0.97276259, 1.2769792 ,
#        1.        , 1.        , 0.97116494, 0.99074606, 1.        ,
#        1.        , 1.        , 1.37019302, 0.9427776 , 0.9427776 ,
#        1.        , 1.        , 1.        , 1.        , 1.        ,
#        0.93384694, 1.02599257, 1.        , 0.23374024, 0.99074606,
#        1.        , 1.01508804, 1.34410843, 0.9427776 , 1.19816956,
#        1.        , 1.01621622, 1.        , 0.97116494, 1.        ,
#        1.        , 0.93384694, 0.93384694, 1.30687331, 1.        ,
#        1.30687331, 1.39128615, 0.93983279, 1.        , 1.        ,
#        1.        , 1.01508804, 1.        , 1.        , 1.        ,
#        1.        , 1.        , 1.        , 1.        , 0.99074606,
#        1.00977676, 1.        , 1.        , 0.99074606, 1.01621622,
#        0.99074606, 1.        , 1.        , 0.93983279, 1.        ,
#        1.        , 1.        , 1.        , 1.        , 0.93384694,
#        1.        , 1.2769792 , 0.97276259, 1.        , 0.99074606,
#        1.        , 0.99074606, 1.        , 1.        , 1.        ,
#        1.2769792 , 1.00977676, 0.9427776 , 1.        , 1.        ,
#        1.        , 1.        , 0.74621224, 1.        , 1.        ,
#        1.02704152, 1.        , 0.99074606, 1.01621622, 0.95970297])
