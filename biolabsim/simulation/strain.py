
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, TYPE_CHECKING

from pandas import DataFrame
from cobra import Model

from .metabolism import load_local_sbml_model, Help_GeneAnnotator, Help_Expr2Flux, Help_FluxCalculator
from .genome import Help_GenomeGenerator, make_GenomeBckgd
from ..manipulation.genetic import Help_MutActProm
from ..common import Sequence, PromoterSite



class Strain (ABC) :
    """
    Abstract class that stores information of genome and other associated metabolic information.
    """

    # Strain ID or name.
    name : str

    # Path where the metabolic network model is stored.
    model_path : Path

    # CobraPy metabolic network model.
    model : Model

    # Information about the genes.
    genes_df : DataFrame # ['RctID','Expression','Promoter','ORF','Fluxes','Expr2Flux']

    # Metabolic model objective score (TODO: improve this description)
    objective : float

    # Genome written in string using the four ACGT bases.
    genome : Sequence


    def get_genome_size (self) -> int :
        return len(self.genome)


    def get_genome_gc_content (self) -> float :
        return round(
            (self.genome.count('G') + self.genome.count('C'))
            / self.get_genome_size()
        , 2)



class WildtypeStrain (Strain) :
    """
    Strains that are considered wildtype which get their genomic information from an SBML Model.
    """

    def __init__ (
        self, name:str, host_name:str, model_path:str,
        genome_size:int = 500, genome_gc_content:float = 0.6
    ) :

        self.name = name
        self.model_path = model_path
        self.model = load_local_sbml_model(self.model_path)

        self.genes_df = Help_GeneAnnotator(host_name, self.model)
        self.genes_df['Fluxes'], self.objective = Help_FluxCalculator(host_name, self)
        self.genes_df['Expr2Flux'] = Help_Expr2Flux(self.genes_df)

        wt_genome = Help_GenomeGenerator(self.genes_df, genome_size, genome_gc_content)
        self.genome = Sequence(wt_genome)



class MutatedStrain (Strain) :
    """
    Strains that mutated from another reference Strain.
    """

    def __init__ (
        self, name:str, host_name:str, ref_strain:Strain,
        num_enzymes:int = 3, target_site:PromoterSite = '-10', num_mutations:int = 2
    ) :

        self.name = name
        self.model_path = ref_strain.model_path
        self.model = ref_strain.model

        # use data of the reference strain to created mutates gene sequences
        mut_genome, mut_genes_df = Help_MutActProm(
            str(ref_strain.genome), ref_strain.genes_df,
            NumberEnzymes=num_enzymes, Target=target_site, NumberMutations=num_mutations
        )
        self.genome = Sequence(mut_genome)

        self.genes_df = mut_genes_df
        self.genes_df['Fluxes'], self.objective = Help_FluxCalculator(host_name, ref_strain, self)
        self.genes_df['Expr2Flux'] = Help_Expr2Flux(self.genes_df)




class FabricatedStrain (Strain) :
    """
    Strains that contain no real genes and whose genome will be randomly generated using the
    nucleic bases and a given GC-Content. These genomes do not follow biological rules, but are
    useful to generate very small genome sequences to test and observe the other methods.
    """

    def __init__ (
        self, name:str, genome_size:int = 500, genome_gc_content:float = 0.6
    ) :

        self.name = name

        # Generate an empty gene dataframe and a random genome sequence.
        self.genes_df = DataFrame(
            columns=['RctID','Expression','Promoter','ORF','Fluxes','Expr2Flux']
        )
        self.genome = make_GenomeBckgd( genome_size, genome_gc_content )