from Fluxer import Fluxer
from Compact import Compact
from SLIP import SLIP

fluxers = {None: Fluxer,
           'compact': Compact,
           'slip': SLIP,
           }

if __name__ == '__main__':
    print fluxers
