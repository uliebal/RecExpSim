"""
An host is a unit that can take in any number of host modules. It serves as the container
of all modules and allows for obverser-event communication between the modules. In addition, it is
capable of properly copying itself and it's modules. It also provides some facilities for
pseudo-random generation to allow deterministic generation.

In technical terms, the Host uses composition-over-inheritance to define all behaviours that
Modules may add to it. In addition to the usual forwarding methods, the Host is also an
observable that can be used by each module to communicate with other modules.
"""

from __future__ import annotations
from abc import ABC, abstractmethod
from typing import List, Callable, NamedTuple, Optional, Type, Dict

# from .experiment import Experiment
from .events import Event
from .random import Generator



ListenerCallback = Callable[[Event],None]



class ListenerEntry ( NamedTuple ) :
    evtype: Type[Event]
    run: ListenerCallback



class Host :

    # modules: Dict[str,Module]
    # def __init__ ( self, exp: Experiment, module_defs: Optional[List[Type[Module]]], ref: Optional[Host] ) :
    #     self.exp = exp
    #     self.listeners = []
    #     if module_defs is not None : # Create new with requested Module types.
    #         for ModDef in module_defs :
    #             self.modules[ ModDef( host=self ) ]
    #     elif ref is not None : # Duplicate given a reference Host.
    #         for mod in ref.modules :
    #             self.modules[ type(mod) ] = mod.clone(self)
    #     else :
    #         raise Exception("Attempted to initialize an Host without a Module definitions nor a reference Host.")

    exp: 'Experiment'

    listeners: List[ListenerEntry]

    rnd_seed: int



    def __init__ ( self, exp: 'Experiment', seed:Optional[int]=None ) :
        """
        Init is responsible to initialize all modules of an host, may it be new from scratch or
        by using another host as a reference.

        Code from derived classes should be:

            super().__init__( exp=exp, ref=ref )

            self.moduleA = Module( host=self )
            self.moduleB = Module( host=self, depModule=self.moduleA )
        """
        exp.bind_host(self)
        self.exp = exp
        self.listeners = []
        self.rnd_seed = seed



    def clone ( self ) -> Host :
        """ Clones a new Host on the same experiment. """
        raise Exception("Cloning not implemented for host '{}'.".format(type(self).__name__))


    def emit ( self, event: Event ) -> None :
        """ Trigger all listeners for a given Event. """
        for le in self.listeners :
            if type(event) == le.evtype :
                le.run(event) # run the listener with the emitted event.



    def observe ( self, evtype: Type[Event], run: ListenerCallback ) -> None :
        """ Bind a new listener for a type of event. """
        self.listeners.append( ListenerEntry( evtype=evtype, run=run ) )



    def make_generator ( self ) -> Generator :
        """ Construct a random number generator with the same seed stored in the host. """
        return Generator(self.rnd_seed)



class HostException (Exception) :
    """ Exception that is triggered when a Host cannot be created. """
    pass
