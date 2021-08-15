"""
This is a recreation of the RecExpSim Notebook.
It should be roughly similar to the existing notebook while using the new restructured codebase.
"""

# Use biolabsim from parent folder.
import os
import sys
sys.path.append( os.path.abspath(os.path.join('.')) )
sys.path.append( os.path.abspath(os.path.join('..')) )


# 1. Set-up of simulation environment
import os
import numpy as np
import matplotlib.pyplot as plt
from silvio.catalog.recexpsim_original import RecExperiment
print('System ready')


# 2. Lab setup
exp = RecExperiment( seed=100, equipment_investment=1500 )
host = exp.create_host( 'ecol' )
print()
host.print_status()


# 3.1. Experiment set-up
temperatures = [26,30,34] # [22,24,26,28,30]
print()
cult_df = host.sim_growth(temperatures)
print(cult_df.to_string())


# 3.2.2.1 Visual analysis of exponential growth
Time, Biomass = cult_df.iloc[:,0], cult_df.iloc[:,1:]
LnBiomass = np.log(Biomass)
[plt.scatter(Time, X, label=Exp) for Exp,X in LnBiomass.iteritems()]
plt.legend()
plt.show()


# 3.2.2.2 Determine maximum biomass and growth rate
Idx_optT, Linear_optT = 0, 7
MB = np.mean(Biomass.iloc[Linear_optT:,Idx_optT])
GR = np.polyfit(Time.iloc[:Linear_optT],LnBiomass.iloc[:Linear_optT,Idx_optT],1)
print()
print('max biomass: {:.0f}\nmax growth rate: {:.2f}'.format(MB,GR[0]))


# 4.1.1 Promoter choice and cloning
from silvio.extensions.records.gene.crafted_gene import CraftedGene
from Bio.Seq import Seq
Promoter1 = Seq("GCCCAAAAAAAAAGCAAACACGTAAAGGAAAAAATGCACG")
Primer1 = Seq("CGGGTTTTTTTTTCG")
Tm = 48 # melting temperature
NewGene = CraftedGene( name="MyGenA", prom=Promoter1, orf="GGGGGGGGGG" )

print("\nTry cloning with bad melting temperature.")
bad_host, clone_outcome = host.clone_with_recombination( Primer1, gene=NewGene, tm=100 )
print("Clone Outcome: " + clone_outcome)
bad_host.print_status()

print("\nTry cloning with good parameters.")
good_host, clone_outcome = host.clone_with_recombination( Primer1, gene=NewGene, tm=Tm )
print("Clone Outcome: " + clone_outcome)
good_host.print_status()


# 4.1.2 Measurement of the promoter strength
prom_str = good_host.calc_promoter_strength( NewGene )
print("\nPromoter Strength of {} is: {}".format(NewGene.name, prom_str))


# 4.1.3 Measurement of the final vaccine expression rate
# in Make_ProductionExperiment: Clone_ID (string), Opt. Temp (int), Opt. Growth rate (float), Opt. Biomass (int)

bad_prod_outcome = good_host.sim_vaccine_production( gene=NewGene, cult_temp=26, growth_rate=0.91, biomass=49 )
print('Final vaccine production rate: {} {}'.format( bad_prod_outcome.value, bad_prod_outcome.error ))

good_prod_outcome = good_host.sim_vaccine_production( gene=NewGene, cult_temp=26, growth_rate=0.91, biomass=49 )
print('Final vaccine production rate: {} {}'.format( good_prod_outcome.value, good_prod_outcome.error ))

# 4.2.1 Visualization of the results
# DataFile = 'Production_Experiments.csv'
# my_data = np.genfromtxt(DataFile, delimiter=',', skip_header=1).reshape(-1,7)

# GCcont, Express = my_data[:,2], my_data[:,6]
# plt.plot(GCcont,Express, linestyle = '--', marker = 'x', color = 'grey')
# plt.gca().set(xlabel='GC-cont', ylabel='rel. expression', xlim=(.4,.8), ylim=(0,1))
# plt.savefig('RelExpress_Vs_GCcont_allProm.png', format='png')
