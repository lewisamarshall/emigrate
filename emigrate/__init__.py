"""Emigrate: An electrophoresis solver.

Class Frame constructs an electrophoresis system.

Class Solver solves the system.

Class FrameSeries creates an on-disk representation of the solution.

CLI is a command line interface for the solver.
"""
from .Solver import Solver
from .Frame import Frame
from .FrameSeries import FrameSeries
from .CLI import CLI
