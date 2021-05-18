"""
Gene crafted by the user itself.
"""

from __future__ import annotations

from Bio.Seq import Seq

from .gene import Gene



class CraftedGene (Gene) :

    name: str

    orf: Seq

    prom: Seq

    def __init__ ( self, name:str, orf:Seq, prom:Seq ) :
        super().__init__()
        self.name = name
        self.orf = orf
        self.prom = prom

    def clone ( self ) -> CraftedGene :
        return CraftedGene(
            name=self.name,
            orf=self.orf, # Can pass because Seq is immutable.
            prom=self.prom, # Can pass because Seq is immutable.
        )

    def get_name ( self ) -> str :
        return self.name

    def get_orf ( self ) -> Seq :
        return self.orf

    def get_prom ( self ) -> Seq :
        return self.prom

    def make_variant ( self, orf_loc:int, orf_sub:Seq ) -> CraftedGene :
        """ TODO: Maybe its better to clone-then-mutate. """
        new_orf = self.orf[:orf_loc-1] + orf_sub + self.orf[orf_loc:] # orf_loc starts at 1
        return CraftedGene(
            name=self.name,
            orf=new_orf,
            prom=self.prom, # Can pass because Seq is immutable.
        )

