
from typing import Optional, Union, Generic, TypeVar
from dataclasses import dataclass

from pandas import DataFrame, Series



T = TypeVar('Value')



class SimulationException (Exception) :
    """ Base class for all simulation exceptions. """
    pass



@dataclass(frozen=True)
class Outcome (Generic[T]) :
    """ Return from operations with any type of result. """
    value: T
    error: Optional[str] = None # TODO: Think about Outcome dataclasses. I cannot make this optional if extended.

    def has_error ( self ) -> bool :
        return self.error is not None

    def succeeded ( self ) -> bool :
        return not self.has_error()



@dataclass(frozen=True)
class DataOutcome :
    """
    Holds the dataframe of a simulation. Has methods to access whether it worked successfully, to
    print the data or to store it in files.
    """
    value: Union[DataFrame,Series]
    error: Optional[str] = None

    def has_error ( self ) -> bool :
        return self.error is not None

    def succeeded ( self ) -> bool :
        return not self.has_error()

    def display ( self ) -> None :
        if self.value is not None :
            print(self.value)
        else :
            print("Outcome is empty.")

    # TODO Implement this.
    def export ( self, filepath ) -> bool :
        return False
