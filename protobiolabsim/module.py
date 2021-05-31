"""
A Module is an extension of behaviour that can be attached to an Organism.
"""

from __future__ import annotations
from abc import ABC, abstractmethod
from typing import Optional

# from ..base import Organism



class Module (ABC) :

    org: 'Organism'



    def __init__ ( self, org:'Organism', ref:Optional[Module] = None ) :
        """
        Init is resposible to build a new Module, whether it is created from scratch or created
        by using another equally-typed Module as a reference. After building, the Module should
        also register itself on the Organism.

        Arguments are usually in the order of:
            org : The host this module belongs to.
            deps* : Zero or more dependent modules this module is using.
            params* : Zero or more params to build this module.

        Extending Modules have the following structure on their __init__ :

            super().__init__(org,ref)          # Call parent init.
            self.depmodule = depmodule         # Add dependent modules.

            self.param1 = param1               # Set the params from the arguments.

            self.org.bind( 'event_type', self.listen_event_type )   # Register listeners.

        """
        self.org = org



    def clone ( self ) -> Module :
        raise Exception("Cloning not implemented for module '{}'.".format(type(self).__name__))
