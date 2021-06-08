"""
The FermentationModel module contains bioprocess models used to simulate fermentation processes.

Top-most is a general FermentationModel Module. Currently it only holds the MonodModel class, which can be used to
simulate Monod kinetics.In later versions, it could also incorporate other types of model, like nonparametric (
machine learning) models.
"""
from __future__ import annotations
from typing import Literal, Dict
from dataclasses import dataclass

from ...common import OperationOutcome
from ...organism import Organism
from ...module import Module
from protobiolabsim.random import set_seed, pick_uniform


def randomize_cond(seed: int = 100, duration: int = 24) -> Dict[str, float]:
    """
    Randomizes "optimal" process conditions for FermentationModel instance.

    Parameters
    ----------
        seed : int
            Seed for random number generator. Default: 100
        duration: int
            Maximum duration of modeled process [h]. Default: 24 h

    Returns
    -------
        conditions: Dict[str, float]
            Randomized process conditions as Dict
    """
    set_seed(seed)
    conditions = {
        'S0': round(pick_uniform(19, 21), 3),       # Initial substrate concentration [g/L]
        'P0': 0,                                    # Initial product concentration [g/L]
        'X0': round(pick_uniform(0.05, 0.3), 3),    # Initial biomass concentration [g/L]'
        'time': range(1, duration+1),               # Process Duration [h]
        'pH': round(pick_uniform(4, 8), 2),         # optimal pH value (no reference)
        'temp': round(pick_uniform(25, 40), 0)      # optimal growth temperature
        # TODO: Event in case optimal process temperature is used in either Recexpsim.RecOrganism or here
        # 'Temp': self.org.growth.opt_growth_temp
    }
    return conditions


def randomize_params(seed: int = 100) -> Dict[str, float]:
    """
    Randomizes "optimal" kinetic parameters for Monodmodel instance.

    Parameters
    ----------
        seed : int
            Seed for random number generator. Default: 100

    Returns
    -------
        params: Dict[str, float]
            Randomized kinetic parameters as Dict
    """
    set_seed(seed)
    params = {
        'u0':   0,                                  # Initial growth rate [h^-1]
        'umax': round(pick_uniform(0.5, 1.1), 3),   # maximal growth rate [h^-1] (default 0.5 - 1.1)
        'Ks':   round(pick_uniform(2, 8), 3),       # Monod substrate affinity constant (default 2 - 8) [g/L]
        'Yx':   round(pick_uniform(0.4, 0.6), 3),   # Yield coefficient for growth on glucose (default 0.4 - 0.6) [g/g]
        'k1':   round(pick_uniform(0.05, 0.2), 3)  # Production rate of Product (default 0.05 - 0.2) [# h^-1]
    }
    return params


@dataclass(frozen=True)
class FermentationOutcome(OperationOutcome):
    """
    The outcome of a simulated fermentation process.
    Contains a dictionary with lists of calculated values like X (biomass), S (substrate) and P (product).
    """
    results: Dict


class FermentationModel(Module):
    """This is a container for all types of fermentation models, be it parametric (e.g. monod kinetics) or
    nonparametric (e.g. neural network, SVM...). Since most (all?) bioprocess models need initial information about
    the process setup, they can be saved in the conditions attribute.
    """
    # TODO: Design class for params and conditions instead of using Dicts?
    # (Starting) Conditions of fermentation process to be modeled, e.g. Temperature, pH, Substrate concentration(s)...
    conditions: Dict[str, float]
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
    This is a class, which can be used to describe fermentation processes which follow monod kinetics. To achieve
    this, additional information like the mode of operation and kinetic parameters are stored.
    """
    # Mode of operation of modeled fermentation process TODO: Could be extended in future work
    operation_mode: Literal['batch', 'fedbatch', 'continuous']
    # Kinetic parameters depending on type of model and organism, e.g. maximal growth rate µmax
    # TODO: Design class for params and conditions instead of using Dicts?
    params: Dict[str, float]
    # (Starting) Conditions of fermentation process to be modeled, e.g. Temperature, pH, Substrate concentration(s)...
    conditions: Dict[str, float]

    def __init__(self, org: Organism, operation_mode: Literal['batch', 'fedbatch', 'continuous'],
                 params: dict = randomize_params(), conditions: dict = None):
        super().__init__(org)
        self.operation_mode = operation_mode
        self.params = params
        if conditions is not None:
            self.conditions = conditions

    def clone(self) -> MonodModel:
        return MonodModel(
            org=self.org,
            operation_mode=self.operation_mode,
            params=self.params,
            conditions=self.conditions
        )


    def __str__(self) -> str:
        """
        Returns
        -------
        description : str
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
        start_params: dict
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
        else:                                   # TODO: add fedbatch & continuous mode
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
            monod_result: FermentationOutcome
                Result of kinetics (X, S, P, µ, rX, rS, rP) as Lists
        """
        params = self.get_start_values()
        rX, rS, rP = [params['rX0']], [params['rS0']], [params['rP0']]
        time, umax, Ks, = self.conditions['time'], self.params['umax'], self.params['Ks']
        Yx, k1, u = self.params['Yx'], self.params['k1'], [self.params['u0']]
        S, P, X = [self.conditions['S0']], [self.conditions['P0']], [self.conditions['X0']]

        for j in time:
            new_u = umax * S[j - 1] / (Ks + S[j - 1])   # difference quotient of µ
            if new_u >= 0:
                u.append(new_u)
            else:
                u.append(0)

            new_rX = u[j - 1] * X[j - 1]                # difference quotient of Biomass
            if new_rX >= 0:
                rX.append(new_rX)
            else:
                rX.append(0)
            X.append(X[j - 1] + rX[j])                  # New [Biomass]

            new_rS = -(rX[j - 1] / Yx)                  # difference quotient of substrate
            if new_rS <= 0:
                rS.append(new_rS)
            else:
                rS.append(0)
            new_S = S[j - 1] + rS[j]
            if new_S < 0:
                S.append(0)                             # New [Substrate]
            else:
                S.append(new_S)

            new_rP = (k1 * u[j]) * X[j]                 # difference quotient of product
            if new_rP >= 0:
                rP.append(new_rP)
            else:
                rP.append(0)

            P.append(P[j - 1] + rP[j])                  # New [Product]

        results = {
                'X':  X,
                'S':  S,
                'P':  P,
                'u':  u,
                'rX': rX,
                'rS': rS,
                'rP': rP
        }
        monod_result = FermentationOutcome(outcome=True, results=results, message=str(self))

        return monod_result
