"""
The Gene component stores the:
 - ORF and Promoter sequences (enhancers in future?)
 - a name as designation
 - (annotations in future - all that what biopython brings to the table)
"""

from __future__ import annotations
from abc import ABC, abstractmethod

from Bio.Seq import Seq

#from ...experiment import Experiment
from ..interface import Record



class Gene ( Record, ABC ) :

    exp: Experiment

    def __init__ ( self, exp: Experiment ) :
        """ When a gene is created it registers itself in the gene registry. """
        exp.gene_reg.insert(self)
        self.exp = exp

    @abstractmethod
    def clone ( self ) -> Gene :
        pass

    @abstractmethod
    def get_name ( self ) -> str :
        pass

    @abstractmethod
    def get_orf ( self ) -> Seq :
        pass

    @abstractmethod
    def get_prom ( self ) -> Seq :
        pass

    @abstractmethod
    def make_variant ( self, orf_loc:int, orf_sub=Seq ) -> Gene :
        pass