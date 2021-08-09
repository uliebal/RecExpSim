"""
All operations that generate random numbers may come from this globally defined random generator.

This allows the user to provide a seed and generate the same sequences on every run as long as
requests are made in the same order.

The methods it exposes are the same of those from a numpy Generator object.
https://numpy.org/doc/stable/reference/random/generator.html#numpy.random.Generator
"""

from typing import Any, List

from numpy.random import Generator, PCG64



# Initialize the generator with a random seed first.
_randgen = Generator(PCG64())



def set_seed ( new_seed:int ) -> None :
    """
    Sets the seed of the global random generator so that all future random picks follow the same
    pseudo-random sequence.
    """
    global _randgen
    _randgen = Generator(PCG64(new_seed))



def pick_choice ( choices:List[Any], weights:List[float] = None ) -> Any :
    """
    Pick a single value out of a provided list. Can receive a non-normalized weights.
    """
    probs = None
    if weights is not None :
        sum_w = sum(list(weights))
        probs = [ w / sum_w for w in list(weights) ]
    return _randgen.choice( choices, p=probs )


def pick_sample ( choices:List[Any], amount:int ) -> List[Any] :
    """
    Pick a sample from choices.
    """
    return _randgen.choice( choices, size=amount, replace=False )


def pick_integer ( low:int, high:int ) -> int :
    """
    Pick a single integer from a range from low (inclusive) to high (exclusive).
    """
    return _randgen.integers( low, high, endpoint=False )


def pick_uniform ( low:float = 0, high:float = 1 ) -> float :
    """
    Pick a single float from an uniform distribution between low (inclusive) and high (exclusive).
    """
    return _randgen.uniform( low, high )


def pick_normal ( mean:float = 0, sd:float = 0 ) -> float :
    """
    Pick a single number from a normal distribution.
    """
    return _randgen.normal( mean, sd )


def pick_exponential ( beta:float ) -> float :
    """
    Pick a single number from an exponential distribution.
    """
    return _randgen.exponential( beta )
