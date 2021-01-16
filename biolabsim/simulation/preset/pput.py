
from ..host import Host
from ..strain import WildtypeStrain
from ..random import pick_integer
from ..config import METABOLIC_MODEL_DIR


class Pput (Host) :
    """
    Specified "Pseudomonas putida" host.
    """

    def __init__ ( self ) :
        Host.__init__( self,

            name= "Pput",

            # related to the values in Ecol, the values were adjusted according to the ratio of the
            # maximum promoter strengths (0.057/0.04) of the optimal sequences (see expression measurement issue).
            # unit: in gDCW/l, source 1: https://onlinelibrary.wiley.com/doi/pdf/10.1002/bit.25474, source 2: https://link.springer.com/article/10.1385/ABAB:119:1:51
            max_biomass= pick_integer(45,145),

            # Obtained from:  http://bigg.ucsd.edu/models/iJN746
            # TODO: Confirm this is the correct model. There are multiple options.
            strain= WildtypeStrain( name= "WT", host_name= "Pput", model_path= METABOLIC_MODEL_DIR / "iJN746.xml" )

        )
