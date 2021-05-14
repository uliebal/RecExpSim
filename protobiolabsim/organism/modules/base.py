"""
A Module is an extension of behaviour that can be attached to an Organism.
"""

from __future__ import annotations
from abc import ABC, abstractmethod

from ..base import Organism



class Module ( ABC ) :

    org: Organism

    def __init__ ( self, org:Organism ) :
        """ Init should register its own hooks on the organism. """
        self.org = org

    @abstractmethod
    def clone ( self, org:Organism ) -> Module :
        pass