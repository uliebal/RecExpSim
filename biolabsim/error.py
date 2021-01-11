


class HostHasNoStrain (Exception) :
    """Raised when a Host has no Strain but wants to use one."""
    def __init__ ( self ) :
        super().__init__("Host has no Strain but wants to perform an action with a Strain.")
