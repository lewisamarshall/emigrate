from __future__ import absolute_import

from .Equilibrator import Equilibrator
from .Fixed import Fixed
from .VariablepH import VariablepH

equilibrators = {None: Equilibrator,
                 'fixed': Fixed,
                 'pH': VariablepH,
                 }
