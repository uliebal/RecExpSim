"""
GenomeSequence is able to output a sequence given the current gene library.
"""

from __future__ import annotations
from typing import List

from Bio.Seq import Seq

from ..base import Organism
from ..events import InitializeEvent
from .base import Module
from .genome_library import GenomeLibrary



class GenomeSequence ( Module ) :

    org: Organism

    genlib: GenomeLibrary

    def __init__ ( self, org:Organism, genlib:GenomeLibrary ) :
        super().__init__(org)
        self.genlib = genlib
        self.org.bind( 'initialize', self.listen_initialize )


    def clone ( self, org:Organism, genlib:GenomeLibrary  ) -> GenomeSequence :
        """ TODO: Review if genlib should be passed here. """
        new_mod = GenomeSequence( org=org, genlib=genlib )
        return new_mod


    def listen_initialize ( self, event:InitializeEvent ) -> None :
        print("GenomeSequence has initialized.")


    def print_sequence ( self ) -> Seq :
        s = Seq("")
        for gene in self.genlib.genes :
            s += Seq("****") + gene.get_prom() + Seq("**") + gene.get_orf()
        return s + Seq("****")

