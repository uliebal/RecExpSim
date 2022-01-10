
from ..host import Host
from ..strain import FabricatedStrain
from ...random import pick_integer
from ...config import METABOLIC_MODEL_DIR


class FabricatedHost (Host) :
    """
    Host that has no real genes, but uses a randomly generated background DNA with specified size
    and GC-content.
    """

    def __init__ ( self, genome_size=80, gc_content=0.6 ) :
        Host.__init__( self,

            name= "Fabrication",

            max_biomass= pick_integer(30,100),

            # Start with a known WT metabolic model.
            strain= FabricatedStrain(
                name= "FabGen",
                genome_gc_content=gc_content,
                genome_size=genome_size
            )

        )

