"""
An organism that implements all currently existing modules.
"""

from __future__ import annotations
from typing import Optional

from ..utils import coalesce
from ..record.gene.base import Gene
from .base import Organism
from .events import InitializeEvent, InsertGeneEvent, RemoveGeneEvent, PromStrCalcEvent
from .modules.genome_library import GenomeLibrary
from .modules.genome_sequence import GenomeSequence
from .modules.phenotype_size import PhenotypeSize
from .modules.prom_str_calc import PromStrCalc



class CompleteOrganism ( Organism ) :

    genlib: GenomeLibrary
    genseq: GenomeSequence
    phenosize: PhenotypeSize
    promstrcalc: PromStrCalc


    def __init__ (
        self,
        exp: Experiment,
        genlib:Optional[GenomeLibrary] = None,
        genseq:Optional[GenomeSequence] = None,
        phenosize:Optional[PhenotypeSize] = None,
        promstrcalc:Optional[PromStrCalc] = None,
    ) :
        """ TODO: Review cloning of genlib dependency on genseq. """
        super().__init__(exp)

        self.genlib = genlib.clone(self) if genlib is not None else GenomeLibrary( org=self )
        self.genseq = genseq.clone(self,self.genlib) if genseq is not None else GenomeSequence( org=self, genlib=self.genlib )
        self.phenosize = phenosize.clone(self) if phenosize is not None else PhenotypeSize( org=self )
        self.promstrcalc = promstrcalc.clone(self,self.genlib) if promstrcalc is not None else PromStrCalc( org=self, genlib=self.genlib )

        event = InitializeEvent()
        self.emit(event)


    def clone ( self ) -> CompleteOrganism :
        """ TODO: Maybe receive Experiment. And use exp.make_clone(org1). """
        return CompleteOrganism(
            exp=self.exp,
            genlib=self.genlib,
            genseq=self.genseq,
            phenosize=self.phenosize,
            promstrcalc=self.promstrcalc,
        )


    def insert_gene ( self, gene: Gene ) -> None :
        event = InsertGeneEvent(gene)
        self.emit(event)


    def replace_gene ( self, old_gene: Gene, new_gene: Gene ) -> None :
        rem_event = RemoveGeneEvent(old_gene)
        self.emit(rem_event)
        ins_event = InsertGeneEvent(new_gene)
        self.emit(ins_event)


    def calc_promoter_strength ( self, gene_name: str ) -> float :
        # find gene by name in genlib
        selected_gene:Optional[Gene] = None
        for g in self.genlib.genes :
            if g.get_name() == gene_name :
                selected_gene = g # gene was found

        result = 0
        if selected_gene is not None :
            result = self.promstrcalc.calculate(selected_gene)

        event = PromStrCalcEvent(gene=selected_gene, result=result)
        self.emit(event)
        return result

