"""
This is a TESTING module to showcase Organism events.
The number inside this module (example: a size of a cell) will increase with the number of genes.
"""

from __future__ import annotations
from typing import List

from ...record.gene.base import Gene
from ..base import Organism
from ..events import InsertGeneEvent, RemoveGeneEvent
from .base import Module



class PhenotypeSize ( Module ) :

    org: Organism

    size: int

    def __init__ ( self, org:Organism ) :
        super().__init__(org)
        self.size = 0
        self.org.bind( 'insert_gene', self.listen_insert_gene )
        self.org.bind( 'remove_gene', self.listen_remove_gene )


    def clone ( self, org:Organism ) -> PhenotypeSize :
        new_mod = PhenotypeSize( org=org )
        return new_mod


    def listen_insert_gene ( self, event:InsertGeneEvent ) -> None :
        self.size += 1
        print("Incremented size by 1 on PhenotypeSize.")


    def listen_remove_gene ( self, event:RemoveGeneEvent ) -> None :
        self.size -= 1
        print("Decremented size by 1 on PhenotypeSize.")
