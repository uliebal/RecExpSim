"""
Contains various utility methods.
"""

def coalesce ( *arg ) :
    """
    Return the first non-null argument given. Emulates a nullish coalesce operator (PEP 505).
    """
    return next( (a for a in arg if a is not None), None )
