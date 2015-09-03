"""An electrophoresis solver."""
# External package imports
from scipy.integrate import ode
import copy

# Emigrate imports
from .flux_schemes import fluxers
from .equilibration_schemes import equilibrators
from .Sequence import Sequence
from .preconditioner import preconditioner


class Solver(object):

    """A class for performing electromigration calculations.

    Args:
        system, path=False, precondition=False
        flux_mode='compact', equilibrium_mode='pH'

    System should be an Frame system that defines the initial condition.

    Filename should be for an hdf5 file to use as the backend for the result.
    If not filename is specified, the solution will be stored in memory.
    """

    # State Information
    initial_condition = None
    state = None

    # Solver Parameters
    _atol = 1e-12
    _rtol = 1e-6
    _method = 'dopri5'

    # Modules
    flux_mode = None
    equilibrium_mode = None
    _fluxer = None
    _equilibrator = None

    def __init__(self, initial_condition, precondition=False,
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

        # Precondition if requested.
        if precondition:
            preconditioner(self.state, self._fluxer)

        # Create empty solution dictionaries
        self.ion_names = [ion.name for ion in self.initial_condition.ions]

    def _set_equilibrium_mode(self):
        """Import an equilibration object to calculate ion properties."""
        # Find the equilibrator.
        try:
            self._equilibrator = equilibrators[self.equilibrium_mode]
        except:
            error_string = '{} is not an equilibrator.'
            raise RuntimeError(error_string.format(self.equilibrium_mode))

        # Initialize the equilibrator.
        self._equilibrator = self._equilibrator(self.state)
        self._equilibrator.equilibrate()

    def _set_flux_mode(self):
        """Import a flux calculator to calculate ion fluxes."""
        # Find the fluxer.
        try:
            self._fluxer = fluxers[self.flux_mode]
        except:
            error_string = '{} is not an fluxer.'
            raise RuntimeError(error_string.format(self.flux_mode))
        # Initialize the fluxer.
        self._fluxer = self._fluxer(self.state)

    def set_reference_frame(self, frame=None, edge='right'):
        """Set the frame of reference.

        Frame should be an ion. Edge should be right or left.
        """
        # #TODO:20 Change this implementation.
        self._fluxer.frame = frame
        self._fluxer.edge = edge

    def solve(self, path='default.hdf5', interval=1., max_time=10):
        """Solve for a series of time points using an ODE solver."""
        if max_time is None:
            raise RuntimeError('Solving requires a finite maximum time.')

        [frame for frame in self.iterate(path, interval, max_time)]

        return Sequence(path, mode='r')

    def iterate(self, path='default.hdf5', interval=1., max_time=None):
        with Sequence(path, mode='w') as sequence:
            sequence.append(self.state)
            self._initialize_solver()
            while self.solver.successful():
                if self.solver.t >= max_time and max_time is not None:
                    return
                else:
                    self._solve_step(interval, max_time)
                    sequence.append(self.state)
                    yield copy.deepcopy(self.state)
            else:
                message = 'Solver failed at time {}.'
                raise RuntimeError(message.format(self.solver.t))

    def _solve_step(self, interval, max_time):
        new_time = min(self.solver.t + interval, max_time)
        self.solver.integrate(new_time)
        self.state.time = self.solver.t
        self._fluxer.unpack(self.solver.y)

    def _initialize_solver(self):
        self._equilibrator.equilibrate()
        self.solver = solver = ode(self._objective)

        solver.set_integrator(self._method, atol=self._atol, rtol=self._rtol)

        if solver._integrator.supports_solout:
            solver.set_solout(self._solout)
        else:
            raise RuntimeError("Solver doesn't support solout."
                               "Equilibrium can't be computed.")

        solver.set_initial_value(self._fluxer.pack(self.initial_condition))

    def _solout(self, time, packed):
        """Perform actions when a successful solution step is found."""
        self._fluxer.unpack(packed)
        self._equilibrator.equilibrate()

    def _objective(self, time, packed):
        """The objective function of the solver."""
        self._fluxer.unpack(packed)
        self._fluxer.update()
        return self._fluxer.pack()
