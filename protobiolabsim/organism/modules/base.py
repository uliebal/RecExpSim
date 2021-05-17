"""
A Module is an extension of behaviour that can be attached to an Organism.
"""

from __future__ import annotations
from abc import ABC, abstractmethod
from typing import Optional

# from ..base import Organism



class Module :

    org: 'Organism'



    def __init__ ( self, org:'Organism', ref:Optional[Module] = None ) :
        """
        Init is resposible to build a new Module, whether it is created from scratch or created
        by using another equally-typed Module as a reference. After building, the Module should
        also register itself on the Organism.

        Extending Modules have the following structure on their __init__ :

            super().__init__(org,ref)          # Call parent init.
            self.depmodule = depmodule         # Add dependent modules.

            if ref is None :                   # Build from scratch ...
                self.genes = set()
            else :                             # ... or build from reference.
                self.genes = copy(ref.genes)

            self.org.bind( 'event_type', self.listen_event_type )   # Register listeners.

        """
        self.org = org



    # def clone ( self, org:Organism ) -> Module :
    #     return type(self)( self.org, self ) # Create a copy using oneself as reference.
