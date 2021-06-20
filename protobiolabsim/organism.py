"""
An organism is a unit that can take in any number of organism modules. It serves as the container
of all modules and allows for obverser-event communication between the modules. In addition, it is
capable of properly copying itself and it's modules. It also provides some facilities for
pseudo-random generation to allow deterministic generation.

In technical terms, the Organism uses composition-over-inheritance to define all behaviours that
Modules may add to it. In addition to the usual forwarding methods, the Organism is also an
observable that can be used by each module to communicate with other modules.
"""

from __future__ import annotations
from abc import ABC, abstractmethod
from typing import List, Callable, NamedTuple, Optional, Type, Dict

from numpy.random import SeedSequence, Generator, PCG64

# from .experiment import Experiment
from .events import Event
from .random import pick_integer



ListenerCallback = Callable[[Event],None]



class ListenerEntry ( NamedTuple ) :
    evtype: Type[Event]
    run: ListenerCallback



class Organism :

    # modules: Dict[str,Module]
    # def __init__ ( self, exp: Experiment, module_defs: Optional[List[Type[Module]]], ref: Optional[Organism] ) :
    #     self.exp = exp
    #     self.listeners = []
    #     if module_defs is not None : # Create new with requested Module types.
    #         for ModDef in module_defs :
    #             self.modules[ ModDef( org=self ) ]
    #     elif ref is not None : # Duplicate given a reference Organism.
    #         for mod in ref.modules :
    #             self.modules[ type(mod) ] = mod.clone(self)
    #     else :
    #         raise Exception("Attempted to initialize an Organism without a Module definitions nor a reference Organism.")

    exp: 'Experiment'

    listeners: List[ListenerEntry]

    seed_sequence: int



    def __init__ ( self, exp: 'Experiment', seed: Optional[int] = None ) :
        """
        Init is responsible to initialize all modules of an organism, may it be new from scratch or
        by using another organism as a reference.

        Code from derived classes should be:

            super().__init__( exp=exp, ref=ref )

            self.moduleA = Module( org=self )
            self.moduleB = Module( org=self, depModule=self.moduleA )
        """
        exp.bind_organism(self)
        self.exp = exp
        self.listeners = []
        self.seed_sequence = SeedSequence(seed)



    def clone ( self ) -> Organism :
        """ Clones a new Organism on the same experiment. """
        raise Exception("Cloning not implemented for organism '{}'.".format(type(self).__name__))


    def emit ( self, event: Event ) -> None :
        """ Trigger all listeners for a given Event. """
        for le in self.listeners :
            if type(event) == le.evtype :
                le.run(event) # run the listener with the emitted event.



    def observe ( self, evtype: Type[Event], run: ListenerCallback ) -> None :
        """ Bind a new listener for a type of event. """
        self.listeners.append( ListenerEntry( evtype=evtype, run=run ) )



    def make_generator ( self ) -> Generator :
        """ Construct a random number generator with the same seed stored in the organism. """
        return Generator(PCG64( self.seed_sequence ))