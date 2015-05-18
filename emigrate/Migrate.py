"""An electrophoresis solver."""
import numpy as np
import scipy.integrate as integrate
import scipy.interpolate
from collections import OrderedDict
from .data_structure import Electrolyte
import warnings
import flux_schemes
import equilibration_schemes
from .data_structure import Electromigration
# pylint: disable=W0212


class Migrate(object):

    """A class for performing electromigration calculations.

    Args:
        system, filename=False, precondition=False
        flux_mode='compact', equilibrium_mode='pH'

    System should be an Electrolyte system that defines the initial condition.

    Filename should be for an hdf5 file to use as the backend for the result.
    If not filename is specified, the solution will be stored in memory.
    """

    import constants

    # Migrate State
    t = 0.
    x = None
    concentrations = None
    ions = None
    area = None

    # Solver Parameters
    atol = 1e-12
    rtol = 1e-6

    # Modules
    equilibrator = None
    flux_calculator = None

    # Which parameters to monitor
    adaptive_grid = False
    area_variation = False

    # Solutions
    electromigration = None

    def __init__(self, system, filename=False, precondition=False,
                 flux_mode='compact', equilibrium_mode='pH'):
        """Initialize with a system from the constructor class."""

        self.system = system

        # Prepare state
        self.ions = self.system.ions
        self._ion_names = [ion.name for ion in self.ions]
        self.x = np.array(system.nodes)
        self.concentrations = np.array(self.system.concentrations)
        self.area = system.area

        # Set equilibrium mode.
        self.equilibrum_mode = equilibrium_mode
        self._set_equilibrium_mode()

        # Set flux mode
        self.flux_mode = flux_mode
        self._set_flux_mode()
        self.flux_calculator.update_ion_parameters(self.equilibrator)

        # Get information from flux calculator
        self.N = self.flux_calculator.N
        self.adaptive_grid = self.flux_calculator.adaptive_grid
        self.area_variation = self.flux_calculator.area_variation

        # Precondition if requested.
        if precondition:
            self.precondition()

        # Create empty solution dictionaries
        self.electromigration = Electromigration(self._ion_names, filename, mode='w')
        self._write_solution(0, self.x, self.area, self.concentrations)

    def _set_equilibrium_mode(self):
        """Import an equilibration object to calculate ion properties."""
        if self.equilibrum_mode == 'fixed':
            self.equilibrator = equilibration_schemes.Fixed
        elif self.equilibrum_mode == 'pH':
            self.equilibrator = equilibration_schemes.Variable_pH
        else:
            raise RuntimeError('Available equlibibrators are "fixed" and "pH".'
                               )

        self.equilibrator = self.equilibrator(self.ions, self.system.pH,
                                              self.concentrations)
        self.equilibrator.equilibrate(self.concentrations)

    def _set_flux_mode(self):
        """Import a flux calculator to calculate ion fluxes."""
        if self.flux_mode == 'compact':
            self.flux_calculator = flux_schemes.Compact
        elif self.flux_mode == 'compact adaptive':
            self.flux_calculator = flux_schemes.CompactAdaptive
        elif self.flux_mode == 'slip':
            self.flux_calculator = flux_schemes.SLIP
        elif self.flux_mode == 'minmod':
            self.flux_calculator = flux_schemes.MinmodLimited
        else:
            raise RuntimeError
        self.flux_calculator = self.flux_calculator(self.system)
        self.flux_calculator.update_ion_parameters(self.equilibrator)

    def set_reference_frame(self, frame=None, edge='right'):
        """Set the frame of reference.

        Frame should be an ion. Edge should be right or left.
        """
        self.flux_calculator.frame = frame
        self.flux_calculator.edge = edge

    def _decompose_state(self, state):
        """Decompose the state into X and concentrations."""
        self.x = state[:self.N]
        if self.area_variation:
            self.area = state[self.N:self.N*2]
            self.concentrations = state[self.N*2:].reshape(self.concentrations.shape)
        else:
            self.concentrations = state[self.N:].reshape(self.concentrations.shape)

    def _compose_state(self, x, area, concentrations):
        """Compose X and concentrations into a state."""
        x = x.flatten()
        concentrations = concentrations.flatten()
        if self.area_variation:
            area = area.flatten()
            state = np.concatenate((x, area, concentrations))
        else:
            state = np.concatenate((x, concentrations))
        return state

    def _write_solution(self, t, x, area, concentrations):
        """Write the current state to solutions."""
        pH = self.equilibrator.pH
        ionic_strength = self.equilibrator.ionic_strength
        current_electrolyte = \
            Electrolyte(dict(nodes=x, ions=self.ions,
                        concentrations=concentrations,
                        pH=pH, ionic_strength=ionic_strength,
                        voltage=self.flux_calculator.V, current_density=self.flux_calculator.j,
                        area = self.area)
                        )
        self.electromigration.add_electrolyte(t, current_electrolyte)

    def solve(self, tmax, dt=1, method='dopri5'):
        """Solve for a series of time points using an ODE solver."""
        self.solver = solver = integrate.ode(self._objective)

        solver.set_integrator(method, atol=self.atol, rtol=self.rtol)

        if solver._integrator.supports_solout:
            solver.set_solout(self._solout)
        else:
            warnings.warn("Solver doesn't support solout.")

        solver.set_initial_value(self._compose_state(self.x, self.area, self.concentrations))

        while solver.successful() and solver.t < tmax:
            tnew = min(solver.t + dt, tmax)
            solver.integrate(tnew)
            self._decompose_state(solver.y)
            self._write_solution(solver.t, self.x, self.area, self.concentrations)
        self._decompose_state(solver.y)

        if not solver.successful():
            print 'solver failed at time', solver.t

    def precondition(self):
        """Precondition the system to place most grid points at regions of change."""
        # set up the interpolator to get the new parameters
        concentration_interpolator = \
            scipy.interpolate.interp1d(self.x,
                                       self.concentrations,
                                       kind='cubic')
        if self.area_variation:
            area_interpolator = \
                scipy.interpolate.interp1d(self.x,
                                           self.area,
                                           kind='cubic')

        # Update the flux calculator
        self.flux_calculator.update(self.x, self.area, self.concentrations)

        # Get the node cost from the flux calculator
        cost = self.flux_calculator.node_cost()
        cost = self.flux_calculator.differ.smooth(cost)

        # The last cost is poorly calculated, so set it to an intermediate value
        cost[-1] = np.median(cost)

        # get the new grid parameters
        self.x = self._precondition_x(cost)
        self.concentrations = concentration_interpolator(self.x)
        if self.area_variation:
            self.area = area_interpolator(self.x)

        # equilibrate the new system.
        self.equilibrator.cH = self.equilibrator.pH = None
        self.equilibrator.equilibrate(self.concentrations)
        self.flux_calculator.update_ion_parameters(self.equilibrator)

    def _precondition_x(self, cost):
        """Precondition the grid based on a cost."""
        new_x = np.cumsum(1/cost)
        new_x -= new_x[0]
        new_x *= max(self.x)/new_x[-1]
        return new_x

    def _solout(self, t, state):
        """Perform actions when a successful solution step is found."""
        self._decompose_state(state)
        self.equilibrator.equilibrate(self.concentrations)
        self.flux_calculator.update_ion_parameters(self.equilibrator)

    def _objective(self, t, state):
        """The objective function of the solver."""
        # Update local parameters
        self.t = t
        self._decompose_state(state)

        # Update the flux calculator and get the relevant parameters
        self.flux_calculator.update(self.x, self.area, self.concentrations)
        dcdt = self.flux_calculator.dcdt
        dxdt = self.flux_calculator.node_flux
        dadt = self.flux_calculator.area_flux

        # Compose them and return to the solver
        dstatedt = self._compose_state(dxdt, dadt, dcdt)
        return dstatedt
