"""
Experiment is the top-most scope. It includes all organisms and registries and acts as a global
namespace.

TODO: See if this is necessary. Can registries be top-level?
  For now: all Organism, Registry and Records fully depend on an Experiment.
"""

from __future__ import annotations
from typing import Type

from .organism.base import Organism
from .registry.interface import Registry
from .record.gene.base import Gene


class Experiment :

    orgs: list[Organism]

    gene_reg: Registry[Gene]

    def __init__ ( self ) :
        self.orgs = []
        self.gene_reg = Registry[Gene]()

    # def create_organism ( self, org_class:Type[Organism], **kwargs ) -> None :
    #     org = org_class( exp=self, **kwargs )
    #     self.orgs.append( org )
