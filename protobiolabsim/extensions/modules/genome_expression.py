"""
TODO: GenomeExpression needs to be reviewed.
"""

from __future__ import annotations
from copy import copy
from protobiolabsim.experiment import Experiment
from typing import Optional, TypedDict, Tuple, List
from collections import namedtuple
from dataclasses import dataclass

import numpy as np
import pandas as pd

from ...common import Outcome
from ...random import pick_uniform, pick_normal
from ...organism import Organism
from ...module import Module
from .genome_library import GenomeLibrary
from ..records.gene.gene import Gene
from ..utils import Help_PromoterStrength
#from .growth_record import Growth




class GenomeExpression ( Module ) :

    # Dependent module Genome Library holding the genes to express.
    genlib: GenomeLibrary

    # Optimal primer length.
    opt_primer_len: int

    # Factor which influences the range of the promoter strength. TODO: These characteristics maybe go elsewhere.
    infl_prom_streng : float
    species_prom_streng : float

    # Path to parameter files. TODO: use Paths.
    regressor_file:str
    addparams_file:str



    def __init__ (
        self, org:Organism, genlib:GenomeLibrary,
        opt_primer_len:int, infl_prom_streng:float, species_prom_streng:float, regressor_file:str, addparams_file:str,
    ) :
        super().__init__(org)

        self.genlib = genlib

        self.opt_primer_len = opt_primer_len
        self.infl_prom_streng = infl_prom_streng
        self.species_prom_streng = species_prom_streng
        self.regressor_file = regressor_file
        self.addparams_file = addparams_file

        # self.org.observe( InsertGeneEvent, self.listen_insert_gene )
        # self.org.observe( RemoveGeneEvent, self.listen_remove_gene )



    def clone ( self, org:Organism ) -> GenomeExpression :
        return GenomeExpression(
            org= org,
            genlib= self.genlib,
            opt_primer_len= self.opt_primer_len,
            infl_prom_streng= self.infl_prom_streng,
            species_prom_streng= self.species_prom_streng,
            regressor_file= self.regressor_file,
            addparams_file= self.addparams_file,
        )



    def calc_prom_str ( self, gene:Gene, ref_prom:str ) -> float :
        final_prom_str = 0

        if gene in self.genlib.genes :
            prom_str = Help_PromoterStrength(
                PromSequence= gene.prom,
                RefPromoter= ref_prom,
                Scaler= 1,
                Similarity_Thresh= .4,
                Regressor_File= self.regressor_file,
                AddParams_File= self.addparams_file,
            )
            final_prom_str = round(prom_str * self.infl_prom_streng, 2)

        return final_prom_str