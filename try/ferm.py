"""
This script uses the methods and classes defined in fermentation_model.py.
"""
from protobiolabsim.catalog.recexpsim import RecOrganism
from protobiolabsim.experiment import Experiment
from protobiolabsim.extensions.modules.fermentation_model import FermentationModel, MonodModel, randomize_cond, \
    randomize_params

seed = 68                                              # Change seed for different conditions and params
exp = Experiment()
org = RecOrganism(exp)
general_model = FermentationModel(org=org)
monod = MonodModel(org=org, operation_mode='batch', conditions=randomize_cond(seed=seed, duration=20),
                   params=randomize_params(seed=seed))
# alternatively
monod2 = MonodModel(org=org, operation_mode='batch')    # default seed is 100, default duration 24

start_params = monod.get_start_values()                 # only used internally
result = monod.calculate_monod()

# TODO: Clone doesnt contain calculated results. Desired behaviour or should results be saved in model instance?
monod3 = monod.clone()                                  # Doesnt contain result!
