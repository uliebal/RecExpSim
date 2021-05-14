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



class PromStrCalc ( Module ) :

    org: Organism

    genlib: GenomeLibrary

    def __init__ ( self, org:Organism, genlib:GenomeLibrary ) :
        super().__init__(org)
        self.genlib = genlib


    def clone ( self, org:Organism, genlib:GenomeLibrary  ) -> PromStrCalc :
        """ TODO: Review if genlib should be passed here. """
        new_mod = PromStrCalc( org=org, genlib=genlib )
        return new_mod


    def calculate ( self, gene: Gene ) -> float :
        prom:Seq = gene.get_prom()
        result = prom.count("C") * 3 + prom.count("G") * 2
        return result

