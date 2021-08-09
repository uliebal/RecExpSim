"""
Experiment is the top-most scope. It includes all hosts and registries and acts as a global
namespace.

TODO: See if this is necessary. Can registries be top-level?
    For now: all Host, Registry and Records fully depend on an Experiment.
"""

from __future__ import annotations
from typing import Set
from abc import ABC

from .host import Host
from .registry import Registry
from .random import Generator



class Experiment (ABC) :

    hosts: Set[Host]

    # The experiment will hold an internal random generator to allow repeatability.
    rnd_gen: Generator



    def __init__ ( self, seed:Optional[Union[int,SeedSequence]]=None ) :
        self.hosts = set()
        self.rnd_gen = Generator(seed)



    def bind_host ( self, host:Host ) -> None :
        self.hosts.add(host)



    # def create_host ( self, host_class:Type[Host], **kwargs ) -> None :
    #     host = host_class( exp=self, **kwargs )
    #     self.hosts.append( host )
