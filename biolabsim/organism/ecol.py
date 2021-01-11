
from ..host import Host
from ..strain import WildtypeStrain
from ..random import pick_integer
from ..config import MODEL_DIR


class Ecol (Host) :
    """
    Specified "Escherichia coli" host.
    """

    def __init__ ( self ) :
        Host.__init__( self,

            name= "Ecol",

            # the limits for Ecol were set as shown below
            # unit: in gDCW/l, source (german): https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=2&cad=rja&uact=8&ved=2ahUKEwjzt_aJ9pzpAhWGiqQKHb1jC6MQFjABegQIAhAB&url=https%3A%2F%2Fwww.repo.uni-hannover.de%2Fbitstream%2Fhandle%2F123456789%2F3512%2FDissertation.pdf%3Fsequence%3D1&usg=AOvVaw2XfGH11P9gK2F2B63mY4IM
            max_biomass= pick_integer(30,100),

            # Start with a known WT metabolic model.
            strain= WildtypeStrain( name= "WT", host_name= "Ecol", model_path= MODEL_DIR / "e_coli_core.xml" )

        )

