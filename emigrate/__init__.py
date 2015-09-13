"""Emigrate: An electrophoresis solver.

Class Frame constructs an electrophoresis system.

Class Solver solves the system.

Class Sequence creates an on-disk representation of the solution.

Emigrate can be used from the command line.
"""

from .Solver import Solver
from .Frame import Frame
from .Sequence import Sequence
from .__version__ import __version__
