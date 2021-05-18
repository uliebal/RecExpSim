"""
An experiment that uses organisms which implements all currently existing modules.
"""

from __future__ import annotations
from typing import Optional, TypedDict, Final

from ..utils import coalesce
from ..experiment import Experiment
from ..registry import Registry
from ..organism import Organism
from ..extensions.genome.events import InsertGeneEvent, RemoveGeneEvent
from ..extensions.genome import GenomeLibrary, GenomeSequence

from ..extensions.genome.records.gene import Gene
from ..extensions.genome.modules.genome_library import GenomeLibrary
from ..extensions.genome.modules.genome_sequence import GenomeSequence
from ..extensions.phenotype.phenotype_size import PhenotypeSize



class CompleteExperiment (Experiment) :

    # TODO: Are Registries really necessary or does garbage collection handle everything?
    #gene_reg: Registry

    def __init__ ( self ) :
        super().__init__()



class CompleteOrganism (Organism) :

    genlib: GenomeLibrary

    genseq: GenomeSequence

    phenosize: PhenotypeSize



    def __init__ ( self, exp: Experiment, ref: Optional[CompleteOrganism] = None ) :
        super().__init__( exp=exp, ref=ref )

        if ref is not None :
            self.genlib = GenomeLibrary( org=self, ref=ref.genlib )
            self.genseq = GenomeSequence( org=self, genlib=self.genlib, ref=ref.genseq )
            self.phenosize = PhenotypeSize( org=self, ref=ref.phenosize )
        else :
            self.genlib = GenomeLibrary( org=self )
            self.genseq = GenomeSequence( org=self, genlib=self.genlib )
            self.phenosize = PhenotypeSize( org=self )



    def insert_gene ( self, gene: Gene ) -> None :
        self.emit( InsertGeneEvent(gene) )



    def replace_gene ( self, old_gene: Gene, new_gene: Gene ) -> None :
        self.emit( RemoveGeneEvent(old_gene) )
        self.emit( InsertGeneEvent(new_gene) )



    # # TODO: Integrate the prom str module
    # def calc_promoter_strength ( self, gene_name: str ) -> float :
    #     # find gene by name in genlib
    #     selected_gene:Optional[Gene] = None
    #     for g in self.genlib.genes :
    #         if g.get_name() == gene_name :
    #             selected_gene = g # gene was found

    #     result = 0
    #     if selected_gene is not None :
    #         result = self.promstrcalc.calculate(selected_gene)

    #     event = PromStrCalcEvent(gene=selected_gene, result=result)
    #     self.emit(event)
    #     return result

