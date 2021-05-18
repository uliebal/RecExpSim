"""
Experiment is the top-most scope. It includes all organisms and registries and acts as a global
namespace.

TODO: See if this is necessary. Can registries be top-level?
    For now: all Organism, Registry and Records fully depend on an Experiment.
"""

from __future__ import annotations
from typing import Set
from abc import ABC

from .organism import Organism
from .registry import Registry



class Experiment (ABC) :

    orgs: Set[Organism]



    def __init__ ( self ) :
        self.orgs = set()



    def bind_organism ( self, org:Organism ) -> None :
        self.orgs.add(org)



    # def create_organism ( self, org_class:Type[Organism], **kwargs ) -> None :
    #     org = org_class( exp=self, **kwargs )
    #     self.orgs.append( org )
