from Differentiate import Differentiate
import numpy as np
import warnings


class _Flux_Base(object):

    """A base class for the flux modules used by emigrate."""

    use_adaptive_grid = False
    boundary_mode = 'characteristic'
    differentiation_method = '6th-Order'
    j = 0.
    E = None
    V = 0
    x = None
    u = 0.
    area = None
    _area = 1.
    concentrations = None
    smoother = False
    nonnegative = True
    frame = None
    edge = 'right'
    pH = None
    I = None
    water_conductivity = None
    water_diffusive_conductivity = None
    mode = 'voltage'

    def __init__(self, system):
        """Initialize the compact flux solver."""

        # Prepare the grid points from the system nodes
        self._prep_domain(system.nodes)

        # Prepare the voltage/current and bulk flow mode from the system.
        self.V = system.voltage
        self.current_density = system.current_density
        if self.V:
            self.mode = 'voltage'
            if self.current_density:
                warnings.warn('System has both current and voltage. Using voltage.')
        else:
            self.mode = 'current'

        self.u = system.u

        # use system area if it exists, otherwise default to _area
        self.area = system.area
        self._area = self.area or self._area

        # Create the differentiation system.
        self.differ = Differentiate(self.N, self.dz,
                                    method=self.differentiation_method,
                                    smoother=self.smoother)

    def _prep_domain(self, nodes):
        self.x = np.array(nodes)
        self.z = np.linspace(min(self.x), max(self.x), len(self.x))
        self.N = self.x.size
        self.dz = self.z[1]-self.z[0]

    def first_derivative(self, x_input):
        """Calculate the first derivative with respect to z."""
        return self.differ.first_derivative(x_input.T).T

    def second_derivative(self, x_input):
        """Calculate the second derivative with respect to z."""
        return self.differ.second_derivative(x_input.T).T

    def set_boundary(self, flux):
        """Set the boundary condition at the domain edges."""
        if self.boundary_mode == 'fixed':
            flux[:, 0] *= 0
            flux[:, -1] *= 0
        elif self.boundary_mode == 'characteristic':
            pass
        return flux

    def dcdt(self, x, concentrations):
        self.x = x
        self.concentrations = concentrations
        dcdt = self._dcdt()
        if self.nonnegative is True:
            dcdt = self.impose_nonnegativity(concentrations, dcdt)
        if self.frame is not None:
            self._update_reference_frame()
        return dcdt

    def impose_nonnegativity(self, concentrations, dcdt):
        dcdt = np.where(np.greater(concentrations, 0), dcdt, np.maximum(0, dcdt))
        return dcdt

    def update_ion_parameters(self, equlibrator):
        self.mobility = equlibrator.mobility
        self.diffusivity = equlibrator.diffusivity
        self.molar_conductivity = equlibrator.molar_conductivity
        self.pH = equlibrator.pH
        self.water_conductivity = equlibrator.water_conductivity
        self.water_diffusive_conductivity = equlibrator.water_diffusive_conductivity

    def _update_reference_frame(self):
        if self.edge == 'right':
            E = self.E[-1]
            pH = self.pH[-1]
        elif self.edge == 'left':
            E = self.E[0]
            pH = self.pH[0]
        else:
            raise RuntimeError('Edge must be left or right.')
        self.u = E * self.frame.effective_mobility(pH)

    def _dcdt(self):
        """Calculate the flux of chemical species."""
        return self.concentrations * 0.
