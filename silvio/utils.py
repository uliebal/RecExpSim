"""
Contains various utility methods.
"""

from typing import Any



def coalesce ( *arg ) -> Any :
    """
    Return the first non-null argument given. Emulates a nullish coalesce operator (PEP 505).
    """
    return next( (a for a in arg if a is not None), None )



def alldef ( *arg ) -> bool :
    """
    Returns true if all arguments are non-None.
    """
    return all([ a is not None for a in arg ])
