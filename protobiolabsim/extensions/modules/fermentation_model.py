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
from ...random import set_seed, pick_uniform
from .growth_behaviour import GrowthBehaviour



OperationMode = Literal['batch', 'fedbatch', 'continuous']

InhibitionTerm = Literal['substrate_inhib', 'product_inhib', 'diauxie', 'diffusion_inhib']

@dataclass
class ConditionSet :
    S0: float
    P0: float
    X0: float
    time: range
    pH: float
    temp: int

@dataclass
class ParamSet :
    u0: float
    umax: float
    Ks: float
    Yx: float
    k1: float

@dataclass
class FermentationResult :
    X: List[float]
    S: List[float]
    P: List[float]
    u: List[float]
    rX: List[float]
    rS: List[float]
    rP: List[float]



def randomize_cond(seed: int = 100, duration: int = 24) -> ConditionSet:
    """
    Randomizes "optimal" process conditions for FermentationModel instance.

    Parameters
    ----------
        seed:
            Seed for random number generator.
        duration:
            Maximum duration of modeled process [h].

    Returns
    -------
        conditions:
            Randomized process conditions
    """
    set_seed(seed)
    conditions = ConditionSet(
        S0=   round(pick_uniform(19, 21), 3),     # Initial substrate concentration [g/L]
        P0=   0,                                  # Initial product concentration [g/L]
        X0=   round(pick_uniform(0.05, 0.3), 3),  # Initial biomass concentration [g/L]'
        time= range(1, duration+1),               # Process Duration [h]
        pH=   round(pick_uniform(3, 7), 1),       # optimal pH value (no reference)
        temp= int(pick_uniform(25, 40))           # optimal growth temperature
        # TODO: Event in case optimal process temperature is added in Recexpsim.RecOrganism -> update Temp here?
        # 'Temp': self.org.growth.opt_growth_temp
    )
    return conditions



def randomize_params(seed: int = 100) -> ParamSet:
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
    params = ParamSet(
        u0=   0,                                  # Initial growth rate [h^-1]
        umax= round(pick_uniform(0.5, 1.1), 3),   # maximal growth rate [h^-1] (default 0.5 - 1.1)
        Ks=   round(pick_uniform(2, 8), 3),       # Monod substrate affinity constant (default 2 - 8) [g/L]
        Yx=   round(pick_uniform(0.4, 0.6), 3),   # Yield coefficient for growth on glucose (default 0.4 - 0.6) [g/g]
        k1=   round(pick_uniform(0.05, 0.2), 3)   # Production rate of Product (default 0.05 - 0.2) [# h^-1]
    )
    return params



class FermentationModel(Module):
    """This is a container for all types of fermentation models, be it parametric (e.g. monod kinetics) or
    nonparametric (e.g. neural network, SVM...). Since most (all?) bioprocess models need initial information about
    the process setup, those can be saved in the conditions attribute.
    """
    # (Starting) Conditions of fermentation process to be modeled, e.g. Temperature, pH, Substrate concentration(s)...
    conditions: ConditionSet

    # Organism, to which the model is attached. NOTE: Existence is implied by Module.
    # org: Organism

    def __init__(self, org: Organism, conditions: ConditionSet):
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

    # Keep a relationship to the Growth module.
    growth: GrowthBehaviour

    # Mode of operation of modeled fermentation process
    # FIXME: (FIXMEs = personal notes) Add fedbatch and continuous mode calculation, possibly further modes
    operation_mode: OperationMode

    # Optional inhibition terms to model different metabolic effects, like product inhibiton, diauxic growth etc.
    inhibition_terms: List[InhibitionTerm]
        # NOTE: Optional Lists should usually be replaced by empty lists.

    # Kinetic parameters depending on type of model and organism, e.g. maximal growth rate µmax
    params: ParamSet

    def __init__(
            self, org: Organism, growth: GrowthBehaviour,
            params: ParamSet,
            conditions: ConditionSet,
            operation_mode: OperationMode = 'batch',
            inhibition_terms: List[InhibitionTerm] = [] ):
        super().__init__( org=org, conditions=conditions )
        self.growth = growth
        self.params = params
        self.operation_mode = operation_mode
        self.inhibition_terms = inhibition_terms

    def clone(self, growth: GrowthBehaviour) -> MonodModel:
        """
        FIXME: Avoid pass-by-reference and copy all values.
        """
        return MonodModel(
            org=self.org,
            growth=growth, # Need to use the cloned Growth module.
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
        u0, Yx, k1 = self.params.u0, self.params.Yx, self.params.k1
        X0, S0 = self.conditions.X0, self.conditions.S0

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

    def calculate_monod(self) -> FermentationResult:
        """
        Calculates monod kinetics for current model instance.

        - X: Cell dry weight (g/L)
        - S: Substrate concentration (g/L)
        - P: Product concentration (g/L)
        - µ: specific growth rate (h^-1)
        - rX, rS, rP: difference quotients (discrete derivative) of X,S and P respectively

        Returns
        -------
            Result of kinetics (X, S, P, µ, rX, rS, rP) as Lists
        """
        start_values = self.get_start_values() # NOTE: If pedantic, this could be a Dataclass.
        rX, rS, rP = [start_values['rX0']], [start_values['rS0']], [start_values['rP0']]

        Ks, Yx, k1, u, umax = self.params.Ks, self.params.Yx, self.params.k1,\
            [self.params.u0], self.params.umax

        time, S, P, X, temp, pH = self.conditions.time, [self.conditions.S0], [self.conditions.P0],\
            [self.conditions.X0], self.conditions.temp, self.conditions.pH

        # Adapt umax/k1 for suboptimal pH/Temp.      FIXME: This isnt scientifically backed (literature search needed)!
        # The greater the distance of temperature from optimum (30), the smaller µmax (biomass growth)
        opt_temp = self.growth.opt_growth_temp
        umax = umax * (1 - (abs(opt_temp - temp) / 100))  # TODO: Get optimal growth rate from organism?
        # The greater the distance of pH from optimum (5), the smaller k1 (product formation rate)
        opt_pH = self.growth.opt_growth_ph
        k1 = k1 * (1 - (abs(opt_pH - pH) / 100))

        # Calculate discrete monod kinetics via difference quotients FIXME: add source
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

        return FermentationResult( X=X, S=S, P=P, u=u, rX=rX, rS=rS, rP=rP )
