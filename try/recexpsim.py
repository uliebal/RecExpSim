# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %% [markdown]
# # Reconstrured RecExpSim Notebook
#
# This is a recreation of the RecExpSim Notebook.
# It should be roughly similar to the existing notebook while using the new restructured codebase.

# %%
# Use biolabsim from parent folder.
import os
import sys
sys.path.append( os.path.abspath(os.path.join('.')) )
sys.path.append( os.path.abspath(os.path.join('..')) )


# %%
# 1. Set-up of simulation environment
import os
import numpy as np
import matplotlib.pyplot as plt

from silvio.outcome import combine_data
from silvio.catalog.recexpsim_original import RecExperiment
print('System ready')


# %%
# 2. Lab setup
exp = RecExperiment( equipment_investment=1500, seed=1 )
host = exp.create_host('ecol')
print()
host.print_status()


# %%
# 3.1. Experiment set-up
temperatures = [22,30,36,40] #[22,24,26,28,30]
print()
growth_out = host.simulate_growth(temperatures)
print(growth_out)


# %%
# 3.2.2.1 Visual analysis of exponential growth
growth_out.display_plot()


# %%
# 3.2.2.2 Determine maximum biomass and growth rate
Idx_optT, Linear_optT = 2, 6
# TODO: Maybe calculations of MB and GR are hidden inside the Outcome.
#   Maybe its intuitive to draw the chosen point in the plot.
cult_df = growth_out.value
Time, Biomass = cult_df.iloc[:,0], cult_df.iloc[:,1:]
LnBiomass = np.log(Biomass)
MB = np.mean(Biomass.iloc[Linear_optT:,Idx_optT])
GR = np.polyfit(Time.iloc[:Linear_optT],LnBiomass.iloc[:Linear_optT,Idx_optT],1)
print()
print('max biomass: {:.0f}\nmax growth rate: {:.2f}'.format(MB,GR[0]))


# %%
# 4.1.1 Promoter choice and cloning
from silvio.extensions.records.gene.crafted_gene import CraftedGene
from Bio.Seq import Seq
Promoter1 = Seq("GCCCAAAAAAAAAGCAAACACGTAAAGGAAAAAATGCACG")
Primer1 =   Seq("CGGGTTTTTTTTTCG")
Tm =        48 # melting temperature
NewGene = CraftedGene( name="MyGenA", prom=Promoter1, orf="GGGGGGGGGG" )

print("\nTry cloning with bad melting temperature.")
# TODO: Melting temperature is not affecting clone chance.
bad_host, clone_outcome = host.clone_with_recombination( Primer1, gene=NewGene, tm=1000 )
print("Clone Outcome: " + clone_outcome)
bad_host.print_status()

print("\nTry cloning first with good parameters.")
good_host_a, clone_outcome = host.clone_with_recombination( Primer1, gene=NewGene, tm=Tm )
print("Clone Outcome: " + clone_outcome)
good_host_a.print_status()

Promoter2 = Seq("GCCCAAAACCAAAGCAAACACGTAAAGGAAAAAATGCACG")
Primer2 =   Seq("CGGGTTTTGGTTTCG")
Tm2 =       52
NewGene2 = CraftedGene( name="MyGenB", prom=Promoter2, orf="GGGGGGGGGG" )

print("\nTry cloning second with good parameters.")
good_host_b, clone_outcome = host.clone_with_recombination( Primer2, gene=NewGene2, tm=Tm2 )
print("Clone Outcome: " + clone_outcome)
good_host_b.print_status()


# %%
# 4.1.2 Measurement of the promoter strength
prom_out_1 = bad_host.measure_promoter_strength( NewGene )
print("\nPromoter Strength of bad host:\n{}".format(prom_out_1))

prom_out_2 = good_host_a.measure_promoter_strength( NewGene )
print("\nPromoter Strength of good host A:\n{}".format(prom_out_2))

prom_out_3 = good_host_b.measure_promoter_strength( NewGene ) # Use NewGene but HostB only has NewGene2
print("\nPromoter Strength of good host B with wrong gene:\n{}".format(prom_out_3))

prom_out_4 = good_host_b.measure_promoter_strength( NewGene2 )
print("\nPromoter Strength of good host B with correct gene:\n{}".format(prom_out_4))

# Join all measurements.
print("\nAll Promoter Strengths combined:")
all_prom_out = combine_data([ prom_out_1, prom_out_2, prom_out_3, prom_out_4 ])
all_prom_out.display_data()


# %%
# 4.1.3 Measurement of the final vaccine expression rate
# in Make_ProductionExperiment: Clone_ID (string), Opt. Temp (int), Opt. Growth rate (float), Opt. Biomass (int)

vac_out_1 = good_host_a.simulate_vaccine_production( gene=NewGene, cult_temp=26, growth_rate=0.86, biomass=39 )
print('Vaccine production rate A, try 1: (error:{})\n{}'.format( vac_out_1.error, vac_out_1.value ))

vac_out_2 = good_host_a.simulate_vaccine_production( gene=NewGene, cult_temp=36, growth_rate=0.92, biomass=39 )
print('Vaccine production rate A, try 2: (error:{})\n{}'.format( vac_out_2.error, vac_out_2.value ))

vac_out_3 = good_host_b.simulate_vaccine_production( gene=NewGene2, cult_temp=36, growth_rate=0.92, biomass=70 )
print('Vaccine production rate B, try 1: (error:{})\n{}'.format( vac_out_3.error, vac_out_3.value ))

vac_out_4 = good_host_b.simulate_vaccine_production( gene=NewGene2, cult_temp=36, growth_rate=0.92, biomass=39 )
print('Vaccine production rate B, try 2: (error:{})\n{}'.format( vac_out_4.error, vac_out_4.value ))

# Its possible to make a custom combine_data for the production experiments.
all_vac_out = combine_data([ vac_out_1, vac_out_2, vac_out_3, vac_out_4 ])
all_vac_out.display_data()
all_vac_out.export_data('Production_Experiments.csv')


# %%
# 4.2.1 Visualization of the results
GCcont, Express = all_vac_out.value['Promoter_GC-content'], all_vac_out.value['Expression_Rate']
plt.plot(GCcont,Express, linestyle = '--', marker = 'x', color = 'grey')
plt.gca().set(xlabel='GC-cont', ylabel='rel. expression', xlim=(.3,.8), ylim=(0,1))
plt.savefig('RelExpress_Vs_GCcont_allProm.png', format='png')
plt.show()
