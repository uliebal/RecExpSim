
from typing import Optional, Any
from dataclasses import dataclass

from pandas import DataFrame


class OperationException (Exception) :
    """ Base class for all simulation exceptions. """
    pass


@dataclass(frozen=True)
class Outcome :
    """ Return from operations with any type of result. """
    value: Any
    error: Optional[str] # TODO: Really think about Outcome dataclasses. I cannot make this optional if extended.

    def has_error ( self ) -> bool :
        return self.error is not None

    def succeeded ( self ) -> bool :
        return not self.has_error()


@dataclass(frozen=True)
class DataOutcome :
    """
    Holds the dataframe of a simulation. Has methods to access whether it worked successfully, to
    print the data or to store it in files. """
    value: DataFrame
    error: Optional[str] # TODO: Really think about Outcome dataclasses. I cannot make this optional if extended.

    def has_error ( self ) -> bool :
        return self.error is not None

    def succeeded ( self ) -> bool :
        return not self.has_error()

    def display ( self ) -> None :
        if value is not None :
            print(value)
        else :
            print("Outcome is empty.")

    # TODO Implement this.
    def export ( self, filepath ) -> bool :
        return false
