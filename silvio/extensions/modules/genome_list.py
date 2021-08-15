"""
GenomeList is a module that stores the multiple genes an Host may have.

"""

from __future__ import annotations

from ...host import Host
from ...module import Module
from ..records.gene.gene import Gene
from ..events import InsertGeneEvent, RemoveGeneEvent



class GenomeList ( Module ) :

    genes: set[Gene]


    def __init__ ( self, host:Host, genes: set[Gene] = set() ) :
        super().__init__(host)

        self.genes = genes

        self.host.observe( InsertGeneEvent, self.listen_insert_gene )
        self.host.observe( RemoveGeneEvent, self.listen_remove_gene )



    def listen_insert_gene ( self, event:InsertGeneEvent ) -> None :
        self.insert_gene( event.gene )
        print( "Added gene={} to the GenomeLibrary.".format(event.gene.get_name()) )



    def listen_remove_gene ( self, event:RemoveGeneEvent ) -> None :
        self.genes.remove( event.gene )
        print( "Removed gene={} from the GenomeLibrary.".format(event.gene.get_name()) )



    def insert_gene ( self, gene:Gene ) :
        """ Insert a gene at a specific location in the sequence. """
        self.genes.add(gene)
