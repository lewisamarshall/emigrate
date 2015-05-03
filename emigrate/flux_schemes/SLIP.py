"""Define the SLIP flux solver."""
import numpy as np
from Flux_Base import _Flux_Base
from Flux_Limiter import Flux_Limiter
# pylint: disable = W0232, E1101


class SLIP(_Flux_Base):

    """A compact flux solver with numerical dissipation and adaptive grid."""

    use_adaptive_grid = True
    boundary_mode = 'characteristic'
    differentiation_method = 'dissipative'
    j = 0
    E = None
    V = 0
    x = None
    u = 0
    NI = 3
    Vthermal = .025
    concentrations = None
    pointwave = 1
    smoother = True
    # limiter = Flux_limiter(minmod)

    # Returnable quantities

    def _dcdt(self):
        self.set_derivatives()
        self.set_node_flux()
        dcdt = (self.electromigration_dcdt() +
                self.advection_dcdt() +
                self.diffusion_dcdt()
                )
        return dcdt

    def set_node_flux(self):
        """Calculate the flux of nodes."""
        flux = self.first_derivative(self.node_cost() *
                                     self.xz)
        flux = self.differ.smooth(flux)
        flux *= self.pointwave
        flux[0, ] = flux[-1, ] = 0.
        self.node_flux = flux

    # Components of dcdt
    def electromigration_dcdt(self):
        flux = self.limit(self.concentrations, self.electromigration_flux())
        dcdt = np.diff(flux, 1)/self.dz
        dcdt = np.pad(dcdt, ((0, 0), (2, 2)), 'constant',
                      constant_values=((0, 0), (0, 0))) / self.xz
        return dcdt

    def diffusion_dcdt(self):
        """Calculate flux due to diffusion."""
        cD = self.diffusivity * self.concentrations
        diffusion = (self.second_derivative(cD) -
                     self.first_derivative(cD)*self.xzz/self.xz)/self.xz**2
        return diffusion

    def advection_dcdt(self):
        advection = ((self.node_flux-self.u) *
                     self.cz /
                     self.xz)
        return advection

    def electromigration_flux(self):
        """Calculate flux due to electromigration."""
        self.set_E()
        uc = self.mobility * self.concentrations
        electromigration = uc * self.E
        return electromigration

    # Flux Limitation
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
            self.mobility + self.node_flux

    def limiter(self, c):
        """temporary implimentation of a flux limiter."""
        diff = np.diff(c, 1)
        v = diff[:, 2:]
        w = diff[:, :-2]
        L = 0.5*(np.sign(v)+np.sign(w))*np.min(np.fabs([v, w]), 0)
        # L = 1./2.*(v + w)*(1.-np.fabs((v-w)/(np.fabs(v)-np.fabs(w))))
        return L

    # Comonents of node movement

    def node_cost(self):
        """Calculate the cost function of each node."""
        self.set_Kag()
        deriv = np.fabs(self.cz)
        cost = deriv / np.nanmax(deriv, 1)[:, np.newaxis]
        cost = np.nanmax(cost, 0) + self.Kag
        return cost

    def set_Kag(self):
        """Set the Kag parameter for spacing of low-gradient grid points."""
        self.Kag = ((self.N-self.NI)/self.NI) * self.Vthermal / abs(self.V)

    # Helper Functions
    def set_derivatives(self):
        self.xz = self.first_derivative(self.x)
        self.xzz = self.second_derivative(self.x)
        self.cz = self.first_derivative(self.concentrations)

    def set_current(self):
        """Calculate the current based on a fixed voltage drop."""
        self.j = self.V/sum(self.dz / self.conductivity())

    def set_E(self):
        """Calculate the electric field at each node."""
        self.set_current()
        self.E = -self.j/self.conductivity()

    def conductivity(self):
        """Calculate the conductivty at each location."""
        conductivity = np.sum(self.molar_conductivity
                              * self.concentrations, 0)
        return conductivity
