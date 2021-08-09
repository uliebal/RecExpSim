"""
This is a TESTING module to showcase Host events.
The number inside this module (example: a size of a cell) will increase with the number of genes.
"""

from __future__ import annotations
from typing import Optional

from ...host import Host
from ...module import Module
from ..genome.events import InsertGeneEvent, RemoveGeneEvent



class PhenotypeSize (Module) :

    host: Host

    size: int



    def __init__ ( self, host:Host, ref:Optional[PhenotypeSize] = None ) :
        super().__init__(host)

        if ref is not None :
            self.size = ref.size
        else :
            self.size = 0

        self.host.observe( InsertGeneEvent, self.listen_insert_gene )
        self.host.observe( RemoveGeneEvent, self.listen_remove_gene )



    def listen_insert_gene ( self, event:InsertGeneEvent ) -> None :
        self.size += 1
        print("Incremented size by 1 on PhenotypeSize.")



    def listen_remove_gene ( self, event:RemoveGeneEvent ) -> None :
        self.size -= 1
        print("Decremented size by 1 on PhenotypeSize.")
