"""
An organism is a unit that can take in any number of organism modules. It serves as the container
of all modules and allows for obverser-event communication between the modules. In addition, it is
capable of properly copying itself and it's modules.

In technical terms, the Organism uses composition-over-inheritance to define all behaviours that
Modules may add to it. In addition to the usual forwarding methods, the Organism is also an
observable that can be used by each module to communicate with other modules.
"""

from __future__ import annotations
from abc import ABC, abstractmethod

from typing import List, Callable, NamedTuple

#from ..experiment import Experiment
from .events import EventType, Event



ListenerCallback = Callable[[Event],None]



class ListenerEntry ( NamedTuple ) :
    typ: EventType
    run: ListenerCallback



class Organism ( ABC ) :

    exp: 'Experiment'

    listeners: List[ListenerEntry]


    def __init__ ( self, exp: 'Experiment' ) :
        """ TODO: Deal with experiment circular dependency. """
        self.exp = exp
        self.listeners = []


    def emit ( self, event: Event ) -> None :
        for le in self.listeners :
            if event.typ == le.typ :
                le.run(event) # run the listener with the emitted event.


    def bind ( self, typ: EventType, run: ListenerCallback ) -> None :
        self.listeners.append( ListenerEntry( typ=typ, run=run ) )


    def clone ( self ) -> Organism :
        return type(self)( exp=self.exp, ref=self )

