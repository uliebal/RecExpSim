"""
GenomeLibrary is a module that stores the multiple genes an Organism may have.

However, it only stores the references to those Genes. The actual storage place of the genes is
in `registry.gene` which has special ways to allow for copy-on-write.
"""

from __future__ import annotations
from copy import copy
from typing import Optional, Final

from ....organism import Organism
from ....module import Module
from ..records.gene import Gene
from ..events import InsertGeneEvent, RemoveGeneEvent



class GenomeLibrary ( Module ) :

    genes: set[Gene]



    def __init__ ( self, org:Organism, ref:Optional[GenomeLibrary] = None ) :
        super().__init__(org,ref)

        if ref is not None :
            self.genes = copy(ref.genes)
        else :
            self.genes = set()

        self.org.observe( InsertGeneEvent, self.listen_insert_gene )
        self.org.observe( RemoveGeneEvent, self.listen_remove_gene )



    def listen_insert_gene ( self, event:InsertGeneEvent ) -> None :
        self.genes.add( event.gene )
        print( "Added gene={} to the GenomeLibrary.".format(event.gene.get_name()) )



    def listen_remove_gene ( self, event:RemoveGeneEvent ) -> None :
        self.genes.remove( event.gene )
        print( "Removed gene={} from the GenomeLibrary.".format(event.gene.get_name()) )
