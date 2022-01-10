# Script to test combining multiple dataframes and series.

# Allow import from parent folder.
import os
import sys
sys.path.append( os.path.abspath(os.path.join('.')) )
sys.path.append( os.path.abspath(os.path.join('..')) )

from pandas import Series

from silvio import DataOutcome, combine_data


a1 = DataOutcome( value=Series({ "name":"A1", "first": 45 }), error=None ) # DataOutcome[Series]
a2 = DataOutcome( value=Series({ "name":"A2", "first": 63 }), error=None ) # DataOutcome[Series]
comb_1 = combine_data([a1,a2]) # DataOutcome[DataFrame]
print( comb_1 )

b1 = DataOutcome( value=Series({ "name":"B1", "second": 23.425 }), error=None )
b2 = DataOutcome( value=Series({ "name":"B2", "second": 0.732 }), error=None )
comb_2 = combine_data([comb_1,b1,b2])
print( comb_2 )
