"""
An experiment that uses organisms which implements all currently existing modules.
"""

from __future__ import annotations
from typing import Optional, List, Tuple
import os

import pandas as pd
from Bio.Seq import Seq

from ..random import pick_integer, pick_uniform
from ..utils import Help_Progressbar
from ..experiment import Experiment
from ..registry import Registry
from ..organism import Organism

from ..config import DATADIR
from ..extensions.modules.growth_behaviour import GrowthBehaviour
from ..extensions.modules.genome_list import GenomeList
from ..extensions.modules.genome_expression import GenomeExpression
from ..extensions.records.gene.gene import Gene
from ..extensions.utils import check_primer_integrity_and_recombination



class RecExperiment (Experiment) :

    suc_rate: float

    def __init__ ( self ) :
        super().__init__()
        self.suc_rate = pick_uniform( 0.8, 1.0 ) # ErrorRate(EquipInvest, self._Mutant__Resources)



class RecOrganism (Organism) :

    growth: GrowthBehaviour

    genlist: GenomeList

    genexpr: GenomeExpression


    def __init__ (
        self, exp: Experiment,
        opt_growth_temp = None, infl_prom_streng = None, opt_primer_len = None, max_biomass = None,
        regressor_file = None, addparams_file = None,
        species_prom_streng = None
    ) :
        super().__init__(exp=exp)

        self.genlist = GenomeList( org=self )
        self.genexpr = GenomeExpression( org=self, genlib=self.genlist,
            opt_primer_len=opt_primer_len, infl_prom_streng=infl_prom_streng,
            species_prom_streng=species_prom_streng,
            regressor_file=regressor_file, addparams_file=addparams_file
        )
        self.growth = GrowthBehaviour( org=self,
            genexpr=self.genexpr,
            opt_growth_temp=opt_growth_temp, max_biomass=max_biomass
        )


    def clone ( self ) -> RecOrganism :
        """
        TODO: Not sure how I like this cloning. I need to get into deep properties like
            `self.growth.opt_growth_temp` to re-create the modules. I want to call `module.clone`.
        """
        return RecOrganism(
            exp=self.exp,
            opt_growth_temp= self.growth.opt_growth_temp,
            infl_prom_streng= self.genexpr.infl_prom_streng,
            species_prom_streng= self.genexpr.species_prom_streng,
            opt_primer_len= self.genexpr.opt_primer_len,
            max_biomass= self.growth.max_biomass,
            regressor_file= self.genexpr.regressor_file,
            addparams_file= self.genexpr.addparams_file,
        )



    def print_status ( self ) -> None :
        print("Organism Information:")
        print("  opt_growth_temp = {}".format( self.growth.opt_growth_temp ))
        print("  max_biomass = {}".format( self.growth.max_biomass ))
        print("  opt_primer_len = {}".format( self.genexpr.opt_primer_len ))
        print("  Gene Library: {} genes".format(len(self.genlist.genes)))
        for gene in self.genlist.genes :
            print("  - {} = {} * {}".format(gene.name, gene.prom, gene.orf))



    def sim_growth ( self, temps:List[int] ) -> pd.DataFrame :
        """ Simulate a growth under multiple temperatures and return expected biomasses over time. """
        ( df, pauses ) = self.growth.Make_TempGrowthExp(CultTemps=temps)

        wait = 0.001 # has to be adjusted, waiting time for loading bar
        for pause in pauses :
            loading_time = wait * pause.loading_len
            Help_Progressbar(45, loading_time, pause.exp)

        return df



    def clone_with_recombination ( self, primer:Seq, gene:Gene, tm:int ) -> Tuple[Organism,str] :
        """
        Clone the organism while trying to insert a gene. The Clone might or might not contain
        the added Gene.
        Returns: ( cloned_organism, outcome_message )
        """

        # A clone is always made.
        cloned_org:RecOrganism = self.clone()

        ref_prom = 'GCCCATTGACAAGGCTCTCGCGGCCAGGTATAATTGCACG'
        primer_integrity = check_primer_integrity_and_recombination(gene.prom, primer, tm, ref_prom, self.genexpr.opt_primer_len)
        if not primer_integrity.succeeded :
            return ( cloned_org, "Primer Failed: " + primer_integrity.error )

        # Both primer integrity and recombination succeeded. Insert the new gene into the clone.
        cloned_org.genlist.insert_gene(gene)
        return ( cloned_org, "Cloning with recombination succeeded." )



    def calculate_promoter_strength ( self, gene:Gene ) -> float :
        ref_prom = 'GCCCATTGACAAGGCTCTCGCGGCCAGGTATAATTGCACG'
        return self.genexpr.calc_prom_str( gene, ref_prom )



    def produce_vaccine ( self, gene:Gene, cult_temp:int, growth_rate:float, biomass:int ) -> Tuple[float,str] :
        ref_prom = 'GCCCATTGACAAGGCTCTCGCGGCCAGGTATAATTGCACG' # TODO: Too many times this ref_prom is passed.
        if pick_uniform(0,1) > self.exp.suc_rate :
            outcome = self.growth.Make_ProductionExperiment( gene, cult_temp, growth_rate, biomass, ref_prom, accuracy_Test=.9 )
            if outcome.succeeded() :
                return [ outcome.value, 'Experiment succeeded.' ]
            else :
                return [ 0, outcome.error ]
        else:
            return [ 0, 'Experiment failed, bad equipment.' ]



class Ecol (RecOrganism) :
    def __init__ ( self, exp ) :
        opt_growth_temp = pick_integer(25,40)  # unit: degree celsius, source: https://application.wiley-vch.de/books/sample/3527335153_c01.pdf
        infl_prom_streng = pick_integer(30,50) # explanation see Plot_ExpressionRate
        opt_primer_len = pick_integer(16,28) # unit: nt, source: https://link.springer.com/article/10.1007/s10529-013-1249-8
        max_biomass = pick_integer(30,100) # unit: in gDCW/l, source (german): https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=2&cad=rja&uact=8&ved=2ahUKEwjzt_aJ9pzpAhWGiqQKHb1jC6MQFjABegQIAhAB&url=https%3A%2F%2Fwww.repo.uni-hannover.de%2Fbitstream%2Fhandle%2F123456789%2F3512%2FDissertation.pdf%3Fsequence%3D1&usg=AOvVaw2XfGH11P9gK2F2B63mY4IM
        regressor_file = DATADIR / 'expression_predictor' / 'Ecol-Promoter-predictor.pkl'
        addparams_file = DATADIR / 'expression_predictor' / 'Ecol-Promoter-AddParams.pkl'
        super().__init__(
            exp=exp,
            opt_growth_temp=opt_growth_temp,
            infl_prom_streng=infl_prom_streng,
            opt_primer_len=opt_primer_len,
            max_biomass=max_biomass,
            regressor_file=regressor_file,
            addparams_file=addparams_file,
            species_prom_streng=0.057,
        )



class Pput (RecOrganism) :
    def __init__ ( self, exp ) :
        opt_growth_temp = pick_integer(25,40)  # unit: degree celsius, source: https://application.wiley-vch.de/books/sample/3527335153_c01.pdf
        infl_prom_streng = pick_integer(30,50) # explanation see Plot_ExpressionRate
        opt_primer_len = pick_integer(16,28) # unit: nt, source: https://link.springer.com/article/10.1007/s10529-013-1249-8
        max_biomass = pick_integer(45,145) # unit: in gDCW/l, source 1: https://onlinelibrary.wiley.com/doi/pdf/10.1002/bit.25474, source 2: https://link.springer.com/article/10.1385/ABAB:119:1:51
        regressor_file = DATADIR / 'expression_predictor' / 'Ptai-Promoter-predictor.pkl'
        addparams_file = DATADIR / 'expression_predictor' / 'Ptai-Promoter-AddParams.pkl'
        super().__init__(
            exp=exp,
            opt_growth_temp=opt_growth_temp,
            infl_prom_streng=infl_prom_streng,
            opt_primer_len=opt_primer_len,
            max_biomass=max_biomass,
            regressor_file=regressor_file,
            addparams_file=addparams_file,
            species_prom_streng=0.04
        )
