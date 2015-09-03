"""Define the SLIP flux solver."""
import numpy as np
import warnings

from .Fluxer import Fluxer
from .Flux_Limiter import Flux_Limiter


class SLIP(Fluxer):

    """A compact flux solver with numerical dissipation and adaptive grid."""

    # #TODO:60 hide internal functions
    boundary_mode = 'fixed'
    differentiation_method = 'dissipative'
    area_variation = False

    # Tuning for adaptive grid
    NI = 10
    Vthermal = .025
    pointwave = 1
    # limiter = Flux_limiter(minmod)

    def _update(self):
        self.set_derivatives()
        self.set_E()
        self.set_node_flux()
        self.set_area_flux()
        self.dcdt = (self.electromigration_dcdt() +
                     self.advection_dcdt() +
                     self.diffusion_dcdt()
                     )
        self.dcdt[:, 0] = self.dcdt[:, -1] = 0

    def set_area_flux(self):
        if self.area_variation:
            self.area_flux = (self.node_flux-self._frame_velocity) * self.ax
            # self.area_flux[0] = self.area_flux[-1] = 0.
        else:
            self.area_flux = np.zeros(self.state.nodes.shape)

    def set_node_flux(self):
        """Calculate the flux of nodes."""
        flux = self.first_derivative(self.node_cost() * self.xz)
        flux = self.differ.smooth(flux)
        flux *= self.pointwave
        flux[0, ] = flux[-1, ] = 0.
        self.node_flux = flux

    # Components of dcdt
    def electromigration_dcdt(self):
        flux = self.limit(self.state.concentrations,
                          self.electromigration_flux())
        dcdt = np.diff(flux, 1)/self.dz
        dcdt = np.pad(dcdt, ((0, 0), (2, 2)), 'constant',
                      constant_values=((0, 0), (0, 0))) / self.xz
        return dcdt

    def diffusion_dcdt(self):
        """Calculate flux due to diffusion."""
        cD = self.state.diffusivity * self.state.concentrations
        diffusion = (self.second_derivative(cD) -
                     self.first_derivative(cD) *
                     self.xzz/self.xz)/self.xz**2
        return diffusion

    def advection_dcdt(self):
        advection_speed = (self.node_flux -
                           (self.state.bulk_flow - self._frame_velocity))
        advection = advection_speed * self.cz / self.xz
        return advection

    def electromigration_flux(self):
        """Calculate flux due to electromigration."""
        return (self.state.mobility *
                self.state.concentrations *
                self.state.field)

    # Flux Limitation
    def limit(self, concentrations, flux):
        self.set_alpha()
        limited_flux = 0.5*(flux[:, 3:] + flux[:, :-3]) +\
            self.alpha[:, 1:-2] * (np.diff(concentrations, 1)[:, 1:-1] -
                                   self.limiter(concentrations))
        return limited_flux

    def set_alpha(self):
        """Set alpha for dissipation."""
        self.set_characteristic()
        self.alpha = 0.5 * np.maximum(np.fabs(self.characteristic /
                                      self.xz), 0)

    def set_characteristic(self):
        """Calculate the characteristic speed of paramters."""
        self.characteristic = (self.state.bulk_flow - self._frame_velocity) + \
            self.state.field * self.state.mobility - self.node_flux

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
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            deriv = np.fabs(self.cz)
            cost = deriv / np.nanmax(deriv, 1)[:, np.newaxis]
            cost = np.nanmax(cost, 0) + self.Kag()
        return cost

    def Kag(self):
        """Set the Kag parameter for spacing of low-gradient grid points."""
        return ((self.N-self.NI) / self.NI *
                self.Vthermal / abs(self.state.voltage))

    # Helper Functions
    def set_derivatives(self):
        self.xz = self.first_derivative(self.state.nodes)
        self.xzz = self.second_derivative(self.state.nodes)
        self.cz = self.first_derivative(self.state.concentrations)
        if self._area.size > 1:
            self.az = self.first_derivative(self._area)
            self.ax = self.az/self.xz
        else:
            self.az = self.ax = 0.

    def set_current(self):
        """Calculate the current based on a fixed voltage drop."""
        self.state.current = (self.state.voltage /
                              sum(self.dz / self.conductivity() / self._area))
        self.state.current_density = self.state.current/self._area

    def set_E(self):
        """Calculate the electric field at each node."""
        if self.mode is 'voltage':
            self.set_current()
            self._set_field()

        elif self.mode is 'current':
            self.state.current_density = self.state.current/self._area
            self._set_field()
            self.state.voltage = np.sum((self.state.field[:-1] +
                                         self.state.field[1:]) /
                                        2 * np.diff(self.state.nodes))
        else:
            raise RuntimeError()

    def _set_field(self):
        self.state.field = -(self.state.current_density +
                             self.diffusive_current())/self.conductivity()

    def conductivity(self):
        """Calculate the conductivty at each location."""
        conductivity = np.sum(self.state.molar_conductivity *
                              self.state.concentrations, 0)
        conductivity += self.state.water_conductivity
        return conductivity

    def diffusive_current(self):
        """Calculate the diffusive current at each location."""
        diffusive_current = self.first_derivative(
            np.sum(self.state.molar_conductivity/self.state.mobility *
                   self.state.diffusivity*self.state.concentrations,
                   0) + self.state.water_diffusive_conductivity
            )
        return diffusive_current

    def pack(self, frame=None):
        if frame:
            frame.area = np.array(frame.area)
            if frame.area.size != self.N:
                frame.area = np.resize(frame.area, [self.N])
            queued = (frame.nodes.ravel(),
                      frame.area.ravel(),
                      frame.concentrations.ravel(),
                      )
        else:
            queued = (self.node_flux.ravel(),
                      self.area_flux.ravel(),
                      self.dcdt.ravel(),
                      )
        return np.concatenate(queued)

    def unpack(self, packed):
        self.state.nodes = packed[:self.N]
        self.state.area = packed[self.N:self.N*2]
        self.state.concentrations = \
            packed[self.N*2:].reshape(self.state.concentrations.shape)
