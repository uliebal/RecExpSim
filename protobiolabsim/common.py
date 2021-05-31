
from typing import Optional, Any
from dataclasses import dataclass



@dataclass(frozen=True)
class OperationOutcome :
    """ Return from operations without results. """
    outcome: bool
    message: Optional[str]

