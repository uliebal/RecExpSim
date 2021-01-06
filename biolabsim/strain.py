
from __future__ import annotations
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, TYPE_CHECKING

import pandas

from .metabolism import load_local_sbml_model, Help_GeneAnnotator, Help_Expr2Flux, Help_FluxCalculator
from .genome import Help_GenomeGenerator, Help_MutActProm

if TYPE_CHECKING: # Avoid circular dependencies because of typing.
    from .host import Host
    from .common import PromoterSite



class Strain (ABC) :
    """
    Abstract class that stores information of genome and other associated metabolic information.
    """

    # Strain ID or name.
    name : str

    # Path where the metabolic network model is stored.
    model_path : Path

    # CobraPy metabolic network model.
    # TODO: Add cobra internal model typing
    model : Any

    # Information about the genes.
    genes_df : pandas.DataFrame

    # Genome written in string using the four ACTG bases.
    genome : str


    def get_genome_size (self) -> int :
        return len(self.genome)


    def get_genome_gc_content (self) -> float :
        return round((self.genome.count('G') + self.genome.count('C'))/self.genome_size,2)



class WildtypeStrain (Strain) :
    """
    Strains that are considered wildtype which get their genomic information from an SBML Model.
    TODO: Is "genomic" the correct word?
    TODO: Still needs the Host as argument. Might be strange if Strain is a subtype of Host.
    """

    def __init__ (
        self, name:str, host:Host, model_path:str,
        wanted_genome_size:int = 500, wanted_genome_gc_content:float = 0.6
    ) :

        self.name = name
        self.model_path = model_path
        self.model = load_local_sbml_model(self.model_path)

        self.genes_df = Help_GeneAnnotator(host.name, self.model)
        self.genes_df['Fluxes'], self.objective = Help_FluxCalculator(self, host)
        self.genes_df['Expr2Flux'] = Help_Expr2Flux(self.genes_df)

        self.genome = Help_GenomeGenerator(self.genes_df, wanted_genome_size, wanted_genome_gc_content)



class MutatedStrain (Strain) :
    """
    Strains that mutated from another reference Strain.
    TODO: Still needs the Host as argument. Might be strange if Strain is a subtype of Host.
    """

    def __init__ (
        self, name:str, host:Host, ref_strain:Strain,
        num_enzymes:int = 3, target_site:PromoterSite = '-10', num_mutations:int = 2
    ) :

        self.name = name
        self.model_path = ref_strain.model_path
        self.model = ref_strain.model

        # use data of the reference strain to created mutates gene sequences
        mut_genome, mut_genes_df = Help_MutActProm(
            ref_strain.genome, ref_strain.genes_df,
            NumberEnzymes=num_enzymes, Target=target_site, NumberMutations=num_mutations
        )

        self.genes_df = mut_genes_df
        self.genes_df['Fluxes'], self.objective = Help_FluxCalculator(self, host, True)
        self.genes_df['Expr2Flux'] = Help_Expr2Flux(self.genes_df)

        self.genome = mut_genome

