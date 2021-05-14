"""
GenomeLibrary is a module that stores the multiple genes an Organism may have.

However, it only stores the references to those Genes. The actual storage place of the genes is
in `registry.gene` which has special ways to allow for copy-on-write.
"""

from __future__ import annotations
from copy import copy

from ...record.gene.base import Gene
from ..base import Organism
from ..events import InsertGeneEvent, RemoveGeneEvent
from .base import Module



class GenomeLibrary ( Module ) :

    org: Organism

    genes: set[Gene]

    def __init__ ( self, org:Organism ) :
        super().__init__(org)
        self.genes = set()
        self.org.bind( 'insert_gene', self.listen_insert_gene )
        self.org.bind( 'remove_gene', self.listen_remove_gene )


    def clone ( self, org:Organism ) -> GenomeLibrary :
        new_mod = GenomeLibrary( org=org )
        new_mod.genes = copy(self.genes)
        return new_mod


    def listen_insert_gene ( self, event:InsertGeneEvent ) -> None :
        self.genes.add( event.gene )
        print( "Added gene={} to the GenomeLibrary.".format(event.gene.get_name()) )


    def listen_remove_gene ( self, event:RemoveGeneEvent ) -> None :
        self.genes.remove( event.gene )
        print( "Removed gene={} from the GenomeLibrary.".format(event.gene.get_name()) )
