"""
The FermentationModel module contains bioprocess models used to simulate fermentation processes.

Top-most is a general FermentationModel Module. Currently it only holds the MonodModel class, which can be used to
simulate Monod kinetics.In later versions, it could also incorporate other types of model, like nonparametric (
machine learning) models.
"""
from __future__ import annotations
from typing import List, Literal, Dict, TypedDict, Optional
from dataclasses import dataclass

from ...common import OperationOutcome
from ...organism import Organism
from ...module import Module
from protobiolabsim.random import set_seed, pick_uniform


def randomize_cond(seed: int = 100, duration: int = 24) -> TypedDict('conditions', {'S0': float, 'P0': float,
                                                                                    'X0': float, 'time': range,
                                                                                    'pH': float, 'temp': int}):
    """
    Randomizes "optimal" process conditions for FermentationModel instance.

    Parameters
    ----------
        seed:
            Seed for random number generator. Default: 100
        duration:
            Maximum duration of modeled process [h]. Default: 24 h

    Returns
    -------
        conditions:
            Randomized process conditions
    """
    set_seed(seed)
    conditions = {
        'S0':   round(pick_uniform(19, 21), 3),     # Initial substrate concentration [g/L]
        'P0':   0,                                  # Initial product concentration [g/L]
        'X0':   round(pick_uniform(0.05, 0.3), 3),  # Initial biomass concentration [g/L]'
        'time': range(1, duration+1),               # Process Duration [h]
        'pH':   round(pick_uniform(3, 7), 1),       # optimal pH value (no reference)
        'temp': int(pick_uniform(25, 40))           # optimal growth temperature
        # TODO: Event in case optimal process temperature is added in Recexpsim.RecOrganism -> update Temp here?
        # 'Temp': self.org.growth.opt_growth_temp
    }
    return conditions


def randomize_params(seed: int = 100) -> Dict[str, float]:
    """
    Randomizes "optimal" kinetic parameters for Monodmodel instance.

    Parameters
    ----------
        seed:
            Seed for random number generator. Default: 100

    Returns
    -------
        params:
            Randomized kinetic parameters
    """
    set_seed(seed)
    params = {
        'u0':   0,                                  # Initial growth rate [h^-1]
        'umax': round(pick_uniform(0.5, 1.1), 3),   # maximal growth rate [h^-1] (default 0.5 - 1.1)
        'Ks':   round(pick_uniform(2, 8), 3),       # Monod substrate affinity constant (default 2 - 8) [g/L]
        'Yx':   round(pick_uniform(0.4, 0.6), 3),   # Yield coefficient for growth on glucose (default 0.4 - 0.6) [g/g]
        'k1':   round(pick_uniform(0.05, 0.2), 3)   # Production rate of Product (default 0.05 - 0.2) [# h^-1]
    }
    return params


@dataclass(frozen=True)
class FermentationOutcome(OperationOutcome):
    """
    The outcome of a simulated fermentation process.
    Contains a dictionary with lists of calculated values like X (biomass), S (substrate) and P (product).
    """
    results: Dict[str, List[float]]


class FermentationModel(Module):
    """This is a container for all types of fermentation models, be it parametric (e.g. monod kinetics) or
    nonparametric (e.g. neural network, SVM...). Since most (all?) bioprocess models need initial information about
    the process setup, those can be saved in the conditions attribute.
    """
    # TODO: Design class for conditions instead of using Dicts?
    # (Starting) Conditions of fermentation process to be modeled, e.g. Temperature, pH, Substrate concentration(s)...
    conditions: TypedDict('conditions', {'S0': float, 'P0': float, 'X0': float, 'time': range, 'pH': float,
                                         'temp': int})
    # Organism, to which the model is attached
    org: Organism

    def __init__(self, org: Organism, conditions: Dict = randomize_cond()):
        super().__init__(org)

        self.conditions = conditions

    def clone(self) -> FermentationModel:
        return FermentationModel(
                org=self.org,
                conditions=self.conditions,
        )


class MonodModel(FermentationModel):
    """
    This is a class which can be used to describe fermentation processes that follow monod kinetics. To achieve
    this, additional information like the mode of operation and kinetic parameters are stored.
    """
    # Mode of operation of modeled fermentation process
    # FIXME: (FIXMEs = personal notes) Add fedbatch and continuous mode calculation, possibly further modes
    operation_mode: Literal['batch', 'fedbatch', 'continuous']
    # Optional inhibition terms to model different metabolic effects, like product inhibiton, diauxic growth etc.
    # TODO: Design classes for inhibitions and params instead of using Dict/List?
    inhibition_terms: Optional[List[Literal['substrate_inhib', 'product_inhib', 'diauxie', 'diffusion_inhib']]]
    # Kinetic parameters depending on type of model and organism, e.g. maximal growth rate µmax
    params: Dict[str, float]
    # (Starting) Conditions of fermentation process to be modeled, e.g. Temperature, pH, Substrate concentration(s)...
    conditions: TypedDict('conditions', {'S0': float, 'P0': float, 'X0': float, 'time': range, 'pH': float,
                                         'temp': int})

    def __init__(self, org: Organism, operation_mode: Literal['batch', 'fedbatch', 'continuous'] = 'batch',
                 inhibition_terms: List = None, params: Dict = randomize_params(), conditions: Dict = None):
        super().__init__(org)
        self.operation_mode = operation_mode
        if inhibition_terms is not None:
            self.inhibition_terms = inhibition_terms
        self.params = params
        if conditions is not None:
            self.conditions = conditions

    def clone(self) -> MonodModel:
        return MonodModel(
            org=self.org,
            operation_mode=self.operation_mode,
            inhibition_terms=self.inhibition_terms,
            params=self.params,
            conditions=self.conditions
        )

    def __str__(self) -> str:
        """
        Returns
        -------
        description:
            Brief description of operation mode, parameters and conditions used in the model
        """
        description = f'Fermentation model instance for a {self.operation_mode} microbial production process.\n' \
                      f'Conditions: {self.conditions} \nParameters: {self.params}'
        return description

    def get_start_values(self) -> Dict[str, float]:
        """
        This function calculates and returns the initial numeric derivatives of biomass, substrate and product.
        These are only used for the first step of the approximation in calculate_monod.

        Returns
        -------
        start_params:
            Initial process derivatives for first calculation step
        """
        u0, Yx, k1 = self.params['u0'], self.params['Yx'], self.params['k1']
        X0, S0 = self.conditions['X0'], self.conditions['S0']

        if self.operation_mode == 'batch':
            rX0 = u0 * X0                       # initial rate of change for biomass (X)
            if Yx != 0 and S0 != 0:             # Prevent division by zero
                rS0 = -(rX0 / (Yx / S0))        # initial rate of change for substrate (S)
            else:
                rS0 = 0
            rP0 = (k1 * u0) * X0                # initial rate of change for product (P)

        else:                                   # FIXME: add fedbatch & continuous mode
            rX0 = 0
            rS0 = 0
            rP0 = 0

        start_params = {
                'rX0': rX0,
                'rS0': rS0,
                'rP0': rP0
                        }

        return start_params

    def calculate_monod(self) -> FermentationOutcome:
        """
        Calculates monod kinetics for current model instance.

        - X: Cell dry weight (g/L)
        - S: Substrate concentration (g/L)
        - P: Product concentration (g/L)
        - µ: specific growth rate (h^-1)
        - rX, rS, rP: difference quotients (discrete derivative) of X,S and P respectively

        Returns
        -------
            monod_result:
                Result of kinetics (X, S, P, µ, rX, rS, rP) as Lists
        """
        start_values = self.get_start_values()
        rX, rS, rP = [start_values['rX0']], [start_values['rS0']], [start_values['rP0']]

        Ks, Yx, k1, u, umax = self.params['Ks'], self.params['Yx'], self.params['k1'],\
            [self.params['u0']], self.params['umax']

        time, S, P, X, temp, pH = self.conditions['time'], [self.conditions['S0']], [self.conditions['P0']],\
            [self.conditions['X0']], self.conditions['temp'], self.conditions['pH']

        # Adapt umax/k1 for suboptimal pH/Temp.      FIXME: This isnt scientifically backed (literature search needed)!
        # The greater the distance of temperature from optimum (30), the smaller µmax (biomass growth)
        umax = umax * (1 - (abs(30 - temp) / 100))  # TODO: Get optimal growth rate & temperature from organism?
        # The greater the distance of pH from optimum (5), the smaller k1 (product formation rate)
        k1 = k1 * (1 - (abs(5 - pH) / 100))         # TODO: Get optimal pH from organism?

        # Calculate discrete monod kinetics via difference quotients FIXME: add source
        for j in time:
            new_u = umax * S[j - 1] / (Ks + S[j - 1])   # difference quotient of µ
            u.append(max(0, new_u))

            new_rX = u[j - 1] * X[j - 1]                # difference quotient of Biomass
            rX.append(max(0, new_rX))
            X.append(X[j - 1] + rX[j])                  # New [Biomass]

            new_rS = -(rX[j - 1] / Yx)                  # difference quotient of substrate
            rS.append(max(0, new_rS))
            new_S = S[j - 1] + rS[j]
            S.append(max(0, new_S))

            new_rP = (k1 * u[j]) * X[j]                 # difference quotient of product
            rP.append(max(0, new_rP))
            P.append(P[j - 1] + rP[j])                  # New [Product]

        # Save results as FermentationOutcome
        results = {'X':  X, 'S':  S, 'P':  P, 'u':  u, 'rX': rX, 'rS': rS, 'rP': rP}
        # TODO: What is the purpose of outcome and message?
        monod_result = FermentationOutcome(outcome=True, results=results, message=str(self))

        return monod_result
