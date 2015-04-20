from Differentiate import Differentiate
import numpy as np


class _Flux_Base(object):

    """A base class for the flux modules used by emigrate."""

    use_adaptive_grid = False
    boundary_mode = 'characteristic'
    differentiation_method = '6th-Order'
    j = 0
    E = None
    V = 0
    x = None
    concentrations = None
    smoother = False
    nonnegative = True

    def __init__(self, N, dz, V, z):
        """Initialize the compact flux solver."""
        self.N = N
        self.z = z
        self.dz = dz
        self.V = V
        self.differ = Differentiate(N, dz,
                                    method=self.differentiation_method,
                                    smoother=self.smoother)

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
        return dcdt

    def impose_nonnegativity(self, concentrations, dcdt):
        dcdt = np.where(np.greater(concentrations, 0), dcdt, np.maximum(0, dcdt))
        return dcdt

    def update_ion_parameters(self, equlibrator):
        self.mobility = equlibrator.mobility
        self.diffusivity = equlibrator.diffusivity
        self.molar_conductivity = equlibrator.molar_conductivity

    def _dcdt(self):
        """Calculate the flux of chemical species."""
        return self.concentrations * 0.
