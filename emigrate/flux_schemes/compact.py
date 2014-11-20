"""Define the Compact flux solver."""
from Differentiate import Differentiate
import numpy as np

class Compact(object):

    """A compact flux solver with no numerical dissipation or adaptive gird."""

    use_adaptive_grid = False
    j = 0
    E = None
    V = 0
    x = None
    concentrations = None
    boundary_mode = 'characteristic'

    def __init__(self, N, dz, V,  mobility, diffusivity, molar_conductivity):
        """Initialize the compact flux solver."""
        self.N = N
        self.dz = dz
        self.V = V
        self.mobility = mobility
        self.diffusivity = diffusivity
        self.molar_conductivity = molar_conductivity
        self.differ = Differentiate(N, dz, method='6th-Order')

    def first_derivative(self, x_input):
        """Calculate the first derivative with respect to z."""
        return self.differ.first_derivative(x_input.T).T

    def second_derivative(self, x_input):
        """Calculate the second derivative with respect to z."""
        return self.differ.second_derivative(x_input.T).T

    def flux(self, x, concentrations):
        """Calculate the flux of chemical species."""
        self.x = x
        self.concentrations = concentrations
        self.set_E()
        total_flux = self.diffusive_flux() + \
            self.electromigration_flux()
        total_flux = self.set_boundary(total_flux)
        return total_flux

    def set_current(self):
        """Calculate the current based on a fixed voltage drop."""
        self.j = self.V/sum(self.dz / self.conductivity())

    def set_E(self):
        """Calculate the electric field at each node."""
        self.set_current()
        self.E = -self.j/self.conductivity()

    def diffusive_flux(self):
        """Calculate flux due to diffusion."""
        cD = self.diffusivity * self.concentrations
        diffusion = \
            self.second_derivative(cD)
        return diffusion

    def electromigration_flux(self):
        """Calculate flux due to electromigration."""
        uc = self.mobility * self.concentrations
        electromigration = \
            self.first_derivative(uc * self.E)
        return electromigration

    def set_boundary(self, flux):
        """Set the boundary condition at the domain edges."""
        if self.boundary_mode == 'fixed':
            flux[:, 0] *= 0
            flux[:, -1] *= 0
        elif self.boundary_mode == 'characteristic':
            pass
        return flux

    def conductivity(self):
        """Calculate the conductivty at each location."""
        conductivity = np.sum(np.tile(self.molar_conductivity,
                                      (1, self.N))
                              * self.concentrations, 0)
        return conductivity
