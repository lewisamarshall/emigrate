from Fluxer import Fluxer
from Compact import Compact
from CompactAdaptive import CompactAdaptive
from SLIP import SLIP
from minmod_limited import MinmodLimited

fluxers = {None: Fluxer,
           'compact': Compact,
           'adaptive_grid': CompactAdaptive,
           'slip': SLIP,
           'minmod': MinmodLimited,
           }

if __name__ == '__main__':
    print fluxers
