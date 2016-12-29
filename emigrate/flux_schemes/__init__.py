from __future__ import absolute_import
from .Fluxer import Fluxer
from .Compact import Compact
from .SLIP import SLIP

fluxers = {None: Fluxer,
           'compact': Compact,
           'slip': SLIP,
           }
