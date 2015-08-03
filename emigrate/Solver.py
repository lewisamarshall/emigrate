"""An electrophoresis solver."""
# Numerical imports
import numpy as np
from scipy.integrate import ode
from scipy.interpolate import interp1d

# Other Libraries
import warnings
import copy

# Emigrate imports
from flux_schemes import fluxers
from equilibration_schemes import equilibrators
from FrameSeries import FrameSeries
from Frame import Frame
from preconditioner import preconditioner

# pylint: disable=W0212


class Solver(object):

    """A class for performing electromigration calculations.

    Args:
        system, filename=False, precondition=False
        flux_mode='compact', equilibrium_mode='pH'

    System should be an Frame system that defines the initial condition.

    Filename should be for an hdf5 file to use as the backend for the result.
    If not filename is specified, the solution will be stored in memory.
    """

    # State Information
    initial_condition = None
    state = None
    time = None
    frame_series = None

    # Solver Parameters
    _atol = 1e-12
    _rtol = 1e-6
    _method = 'dopri5'

    # Modules
    flux_mode = None
    equilibrium_mode = None
    fluxer = None
    equilibrator = None

    def __init__(self, initial_condition, filename=False, precondition=False,
                 flux_mode='slip', equilibrium_mode='pH'):
        """Initialize with a system from the constructor class."""

        self.initial_condition = initial_condition
        self.state = copy.deepcopy(initial_condition)

        # Set equilibrium mode.
        self.equilibrium_mode = equilibrium_mode
        self._set_equilibrium_mode()

        # Set flux mode
        self.flux_mode = flux_mode
        self._set_flux_mode()
        self.fluxer.update_ion_parameters(self.equilibrator)

        # Precondition if requested.
        if precondition:
            self.initial_condition = preconditioner(self.initial_condition,
                                                    self.fluxer)

        # Create empty solution dictionaries
        ion_names = [ion.name for ion in self.initial_condition.ions]
        self.frame_series = FrameSeries(ion_names, filename, mode='w')

    def _set_equilibrium_mode(self):
        """Import an equilibration object to calculate ion properties."""
        # Find the equilibrator.
        try:
            self.equilibrator = equilibrators[self.equilibrium_mode]
        except:
            error_string = '{} is not an equilibrator.'
            raise RuntimeError(error_string.format(self.equilibrium_mode))

        # Initialize the equilibrator.
        self.equilibrator = self.equilibrator(self.state)
        self.equilibrator.equilibrate()

    def _set_flux_mode(self):
        """Import a flux calculator to calculate ion fluxes."""
        # Find the fluxer.
        try:
            self.fluxer = fluxers[self.flux_mode]
        except:
            error_string = '{} is not an fluxer.'
            raise RuntimeError(error_string.format(self.flux_mode))
        # Initialize the fluxer.
        self.fluxer = self.fluxer(self.state)

    def set_reference_frame(self, frame=None, edge='right'):
        """Set the frame of reference.

        Frame should be an ion. Edge should be right or left.
        """
        self.fluxer.frame = frame
        self.fluxer.edge = edge

    def _write_solution(self):
        """Write the current state to solutions."""
        self.frame_series.add_frame(self.time, self.state)

    def solve(self, interval=1, max_time=10, method='dopri5'):
        """Solve for a series of time points using an ODE solver."""
        if max_time is None:
            raise RuntimeError('Solving requires a finite maximum time.')

        print "Solving..."
        for i in self.iterate(interval, max_time):
            pass
        print "Solved."
        return self.frame_series

    def iterate(self, interval=1., max_time=None):
        self._initialize_solver()
        while self.solver.successful():
            if self.solver.t >= max_time and max_time is not None:
                return
            else:
                self._solve_step(interval, max_time)
                self.t = self.solver.t
                self.fluxer.unpack(self.solver.y, self.state)
                yield self.state
        else:
            message = 'Solver failed at time {}.'
            raise RuntimeError(message.format(self.solver.t))

    def _solve_step(self, dt, tmax):
        tnew = min(self.solver.t + dt, tmax)
        self.solver.integrate(tnew)
        self.fluxer.unpack(self.solver.y, self.state)
        self._write_solution()

    def _initialize_solver(self):
        self.time = 0
        self._write_solution()

        self.solver = solver = ode(self._objective)

        solver.set_integrator(self._method, atol=self._atol, rtol=self._rtol)

        if solver._integrator.supports_solout:
            solver.set_solout(self._solout)
        else:
            warnings.warn("""Solver doesn't support solout.
                             Equilibrium won't be computed.""")

        solver.set_initial_value(self.fluxer.pack(self.initial_condition))

    def _solout(self, time, packed):
        """Perform actions when a successful solution step is found."""
        self.fluxer.unpack(packed, self.state)
        self.equilibrator.equilibrate()
        self.fluxer.update_ion_parameters(self.equilibrator)

    def _objective(self, time, packed):
        """The objective function of the solver."""
        # Update local parameters
        self.time = time
        self.fluxer.unpack(packed, self.state)

        # Update the flux calculator and get the relevant parameters
        self.fluxer.update()
        return self.fluxer.pack()
