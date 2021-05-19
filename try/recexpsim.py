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
#from tabulate import tabulate
from protobiolabsim.catalog.recexpsim import RecExperiment, RecOrganism
print('System ready')


# 2. Lab setup
exp = RecExperiment() # TODO: Budget System not implemented.
host = RecOrganism(exp=exp) # TODO: Use organism with some genes.
host.print_status()


# 3.1. Experiment set-up
temperatures = [22,26,30,34]
cult_df = host.sim_growth(temperatures)
print(cult_df.to_string())





# 3.2.2.1 Visual analysis of exponential growth
# Time, Biomass = my_data[:,0], my_data[:,1:]
# DataFile = 'Strain_characterization_1.csv'
# LnBiomass = np.log(Biomass)
# [plt.scatter(Time, X, label=Exp) for Exp,X in enumerate(LnBiomass.T)]
# plt.legend()
# my_data = np.genfromtxt(DataFile, delimiter=',', skip_header=1)
