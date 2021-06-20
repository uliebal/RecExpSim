"""
This script uses the methods and classes defined in fermentation_model.py.
"""

# Use biolabsim from parent folder.
import os
import sys
sys.path.append( os.path.abspath(os.path.join('.')) )
sys.path.append( os.path.abspath(os.path.join('..')) )

from protobiolabsim.catalog.fermprosim import FermExperiment, SomeAcidicOrganism


seed = 68                                              # Change seed for different conditions and params
exp = FermExperiment()
org = exp.create_acidic_organism(seed)
# alternatively (but less preferred for the educational workflow)
org2 = SomeAcidicOrganism( exp=exp )    # default seed is 100, default duration 24, batch mode, conditions&params randomized

# Run the model.
result = org.calc_monod_kinetics()
# alternatively (but less preferred for the educational workflow)
result2 = org.monod.calculate_monod()

print("ORGANISM:")
org.print_status()
print("RESULT:\n{}".format(result))
print("MODEL:\n{}".format(org.monod))

# TODO: Clone doesnt contain calculated results. Desired behaviour or should results be saved in model instance?
#org3 = org.clone()                                  # Doesnt contain result!
