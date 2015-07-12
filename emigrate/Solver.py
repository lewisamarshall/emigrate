"""An electrophoresis solver."""
# Numerical imports
import numpy as np
from scipy.integrate import ode
from scipy.interpolate import interp1d

# Warning import
import warnings

# Emigrate imports
from flux_schemes import fluxers
from equilibration_schemes import equilibrators
from FrameSeries import FrameSeries
from Frame import Frame

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

    system = None

    # Solver Parameters
    atol = 1e-12
    rtol = 1e-6

    # Modules
    equilibrator = None
    fluxer = None

    # Solutions
    frames = None

    def __init__(self, system, filename=False, precondition=False,
                 flux_mode='slip', equilibrium_mode='pH'):
        """Initialize with a system from the constructor class."""

        self.system = system

        # Prepare state
        self.ions = self.system.ions
        self.x = np.array(system.nodes)
        self.concentrations = np.array(self.system.concentrations)
        self.area = system.area

        # Set equilibrium mode.
        self.equilibrium_mode = equilibrium_mode
        self._set_equilibrium_mode()

        # Set flux mode
        self.flux_mode = flux_mode
        self._set_flux_mode()
        self.fluxer.update_ion_parameters(self.equilibrator)

        # Get information from flux calculator
        self.N = self.fluxer.N
        self.adaptive_grid = self.fluxer.adaptive_grid
        self.area_variation = self.fluxer.area_variation

        # Precondition if requested.
        if precondition:
            self.precondition()

        # Create empty solution dictionaries
        ion_names = [ion.name for ion in self.ions]
        self.frames = FrameSeries(ion_names,
                                  filename,
                                  mode='w')
        self._write_solution(0, self.x, self.area, self.concentrations)

    def _set_equilibrium_mode(self):
        """Import an equilibration object to calculate ion properties."""
        try:
            self.equilibrator = equilibrators[self.equilibrium_mode]
        except:
            error_string = '{} is not an equilibrator.'
            raise RuntimeError(error_string.format(self.equilibrium_mode))

        self.equilibrator = self.equilibrator(self.ions,
                                              self.concentrations)
        self.equilibrator.equilibrate(self.concentrations)

    def _set_flux_mode(self):
        """Import a flux calculator to calculate ion fluxes."""
        try:
            self.fluxer = fluxers[self.flux_mode]
        except:
            error_string = '{} is not an fluxer.'
            raise RuntimeError(error_string.format(self.flux_mode))
        self.fluxer = self.fluxer(self.system)
        self.fluxer.update_ion_parameters(self.equilibrator)

    def set_reference_frame(self, frame=None, edge='right'):
        """Set the frame of reference.

        Frame should be an ion. Edge should be right or left.
        """
        self.fluxer.frame = frame
        self.fluxer.edge = edge

    def _decompose_state(self, state):
        """Decompose the state into X and concentrations."""
        self.x = state[:self.N]
        if self.area_variation:
            self.area = state[self.N:self.N*2]
            self.concentrations = \
                state[self.N*2:].reshape(self.concentrations.shape)
        else:
            self.concentrations = \
                state[self.N:].reshape(self.concentrations.shape)

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
        current_frame = \
            Frame(dict(nodes=x, ions=self.ions,
                       concentrations=concentrations,
                       pH=pH, ionic_strength=ionic_strength,
                       voltage=self.fluxer.V,
                       current_density=self.fluxer.j,
                       area=self.area)
                  )
        self.frames.add_frame(t, current_frame)

    def solve(self, tmax, dt=1, method='dopri5'):
        """Solve for a series of time points using an ODE solver."""
        self.solver = solver = ode(self._objective)

        solver.set_integrator(method, atol=self.atol, rtol=self.rtol)

        if solver._integrator.supports_solout:
            solver.set_solout(self._solout)
        else:
            warnings.warn("Solver doesn't support solout.")

        solver.set_initial_value(self._compose_state(self.x,
                                                     self.area,
                                                     self.concentrations)
                                 )

        while solver.successful() and solver.t < tmax:
            tnew = min(solver.t + dt, tmax)
            solver.integrate(tnew)
            self._decompose_state(solver.y)
            self._write_solution(solver.t,
                                 self.x,
                                 self.area,
                                 self.concentrations)
        self._decompose_state(solver.y)

        if not solver.successful():
            print 'solver failed at time', solver.t

    def precondition(self):
        """Precondition the system by spacing the grid points."""
        # set up the interpolator to get the new parameters
        concentration_interpolator = interp1d(self.x,
                                              self.concentrations,
                                              kind='cubic')
        if self.area_variation:
            area_interpolator = interp1d(self.x,
                                         self.area,
                                         kind='cubic')

        # Update the flux calculator
        self.fluxer.update(self.x, self.area, self.concentrations)

        # Get the node cost from the flux calculator
        cost = self.fluxer.node_cost()
        cost = self.fluxer.differ.smooth(cost)

        # The last cost is poorly calculated, set it to an intermediate value
        cost[-1] = np.median(cost)

        # get the new grid parameters
        self.x = self._precondition_x(cost)
        self.concentrations = concentration_interpolator(self.x)
        if self.area_variation:
            self.area = area_interpolator(self.x)

        # equilibrate the new system.
        self.equilibrator.cH = self.equilibrator.pH = None
        self.equilibrator.equilibrate(self.concentrations)
        self.fluxer.update_ion_parameters(self.equilibrator)

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
        self.fluxer.update_ion_parameters(self.equilibrator)

    def _objective(self, t, state):
        """The objective function of the solver."""
        # Update local parameters
        self.t = t
        self._decompose_state(state)

        # Update the flux calculator and get the relevant parameters
        self.fluxer.update(self.x, self.area, self.concentrations)
        dcdt = self.fluxer.dcdt
        dxdt = self.fluxer.node_flux
        dadt = self.fluxer.area_flux

        # Compose them and return to the solver
        dstatedt = self._compose_state(dxdt, dadt, dcdt)
        return dstatedt
