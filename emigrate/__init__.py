"""Emigrate: An electrophoresis solver.

Class Frame constructs an electrophoresis system.

Class Solver solves the system.

Class Sequence creates an on-disk representation of the solution.

CLI is a command line interface for the solver.
"""
from .Solver import Solver
from .Frame import Frame
from .Sequence import Sequence
