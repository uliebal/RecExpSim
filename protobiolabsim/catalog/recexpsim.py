"""
An experiment that uses organisms which implements all currently existing modules.
"""

from __future__ import annotations
from typing import Optional, List
from random import randint

import pandas as pd

from ..utils import coalesce
from ..experiment import Experiment
from ..registry import Registry
from ..organism import Organism

from ..extensions.physiology.growth_module import GrowthBehaviour



class RecExperiment (Experiment) :

    def __init__ ( self ) :
        super().__init__()



class RecOrganism (Organism) :

    growth: GrowthBehaviour



    def __init__ ( self, exp: Experiment, ref: Optional[RecOrganism] = None ) :
        super().__init__( exp=exp, ref=ref )

        if ref is not None :
            self.growth = GrowthBehaviour( org=self, ref=ref.growth )
        else :
            opt_temp = randint(25,40)
            ecol_max_biomass = randint(30,100)
            self.growth = GrowthBehaviour( org=self, opt_growth_temp=opt_temp, max_biomass=ecol_max_biomass )



    def print_status ( self ) -> None :
        print("Hidden Information:")
        print("  opt_growth_temp = {}".format( self.growth.opt_growth_temp ))
        print("  max_biomass = {}".format( self.growth.max_biomass ))



    def sim_growth ( self, temps:List[int] ) -> pd.DataFrame :
        ( df, pauses ) = self.growth.grow(CultTemps=temps)

        wait = 0.01 # has to be adjusted, waiting time for loading bar
        for pause in pauses :
            loading_time = wait * pause.loading_len
            Help_Progressbar(45, loading_time, pause.exp)

        return df









def Help_Progressbar(n, loading_time, add):
    '''function for display of a loading bar, n: width of loading bar'''
    import sys
    import time

    loading = '.' * n
    for i in range(n+1):
        # this loop replaces each dot with a hash!
        print('\r%s progress{}: %3d percent'.format(add) % (loading, i*100/n), end='')
        loading = loading[:i] + '#' + loading[i+1:]
        time.sleep(loading_time)
    sys.stdout.write("\n")
