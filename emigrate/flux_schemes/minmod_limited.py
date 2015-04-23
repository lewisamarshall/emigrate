"""Define the SLIP flux solver."""
import numpy as np
from Flux_Base import _Flux_Base
from Flux_Limiter import Flux_Limiter
# pylint: disable = W0232, E1101


class MinmodLimited(_Flux_Base):

    """A compact flux solver with numerical dissipation and adaptive grid."""

    use_adaptive_grid = False
    boundary_mode = 'characteristic'
    differentiation_method = 'dissipative'
    j = 0
    E = None
    V = 0
    x = None
    u = 0
    concentrations = None
    # limiter = Flux_limiter(minmod)

    def _dcdt(self):
        self.set_derivatives()
        flux = self.limit(self.concentrations, self.flux())
        dcdt = np.diff(flux, 1)/self.dz
        dcdt = np.pad(dcdt, ((0, 0), (2, 2)), 'constant', constant_values=((0, 0), (0, 0)))
        return dcdt

    def set_derivatives(self):
        self.xz = self.first_derivative(self.x)
        self.xzz = self.second_derivative(self.x)

    def flux(self):
        """Calculate the flux of chemical species."""
        self.set_E()
        total_flux = (self.diffusive_flux() +
                      self.electromigration_flux() +
                      self.advective_flux())
        # total_flux = self.set_boundary(total_flux)
        return total_flux

    def advective_flux(self):
        return -self.u * self.first_derivative(self.concentrations)

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
            self.first_derivative(cD)
        return diffusion

    def electromigration_flux(self):
        """Calculate flux due to electromigration."""
        uc = self.mobility * self.concentrations
        electromigration = uc * self.E
        return electromigration

    def conductivity(self):
        """Calculate the conductivty at each location."""
        conductivity = np.sum(self.molar_conductivity
                              * self.concentrations, 0)
        return conductivity

    def node_flux(self):
        return np.zeros(self.x.shape)

    def limit(self, concentrations, flux):
        self.set_alpha()
        limited_flux = 0.5*(flux[:, 3:] + flux[:, :-3]) +\
            self.alpha[:, 1:-2] * (np.diff(concentrations, 1)[:, 1:-1]
                                   - self.limiter(concentrations))
        return limited_flux

    def set_alpha(self):
        """Set alpha for dissipation."""
        self.set_characteristic()
        self.alpha = 0.5 * np.maximum(np.fabs(self.characteristic /
                                      self.xz), 0)

    def set_characteristic(self):
        """Calculate the characteristic speed of paramters."""
        self.characteristic = self.u + self.E * \
            self.mobility

    def limiter(self, c):
        """temporary implimentation of a flux limiter."""
        diff = np.diff(c, 1)
        v = diff[:, 2:]
        w = diff[:, :-2]
        L = 0.5*(np.sign(v)+np.sign(w))*np.min(np.fabs(np.min([v, w], 1)))
        # L = 1./2.*(v + w)*(1.-np.fabs((v-w)/(np.fabs(v)-np.fabs(w))))
        return L
