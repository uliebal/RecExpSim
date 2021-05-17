"""
GenomeSequence is able to output a sequence given the current gene library.
"""

from __future__ import annotations
from typing import List, Optional

from Bio.Seq import Seq

from ..base import Organism
from ..events import InitializeEvent
from .base import Module
from .genome_library import GenomeLibrary



class GenomeSequence ( Module ) :

    genlib: GenomeLibrary



    def __init__ ( self, org:Organism, genlib:GenomeLibrary, ref:Optional[Module] = None ) :
        super().__init__(org)
        self.genlib = genlib

        if ref is not None :
            pass
        else :
            pass

        self.org.bind( InitializeEvent, self.listen_initialize )



    def listen_initialize ( self, event:InitializeEvent ) -> None :
        print("GenomeSequence has initialized.")



    def print_sequence ( self ) -> Seq :
        s = Seq("")
        for gene in self.genlib.genes :
            s += Seq("****") + gene.get_prom() + Seq("**") + gene.get_orf()
        return s + Seq("****")

