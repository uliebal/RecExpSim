
from typing import Optional
from dataclasses import dataclass


class OperationException (Exception) :
    """ Base class for all simulation exceptions. """
    pass


@dataclass(frozen=True)
class OperationOutcome :
    """ Return from operations without results. """
    success: bool
    error: Optional[str] # TODO: Really think about Outcome dataclasses. I cannot make this optional if extended.

