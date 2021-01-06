
from ..host import Host
from ..random import pick_integer
from ..config import MODEL_DIR


class Pput (Host) :
    """
    Specified "Pseudomonas putida" host.

    TODO: Is this the name?
    """

    def __init__ ( self ) :
        Host.__init__( self,

            name= "Pput",

            # related to the values in Ecol, the values were adjusted according to the ratio of the
            # maximum promoter strengths (0.057/0.04) of the optimal sequences (see expression measurement issue).
            # unit: in gDCW/l, source 1: https://onlinelibrary.wiley.com/doi/pdf/10.1002/bit.25474, source 2: https://link.springer.com/article/10.1385/ABAB:119:1:51
            max_biomass= pick_integer(45,145),

            # TODO: Does it make sense for a Host to have no metabolism? Pput seems to have been
            #   added because of the BioLabSimFun script which does not care about Strains.
            metabolic_model_path= MODEL_DIR / "e_coli_core.xml"

        )



