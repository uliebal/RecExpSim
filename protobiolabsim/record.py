
from __future__ import annotations
from abc import ABC, abstractmethod


class Record ( ABC ) :
    """ Record is a piece of information that can be indexed in a Registry. """

    @abstractmethod
    def clone ( self ) -> Record :
        pass

