"""Define the Compact flux solver."""
import numpy as np
from Fluxer import Fluxer
# pylint: disable = W0232, E1101


class Compact(Fluxer):

    """A compact flux solver with no numerical dissipation or adaptive gird."""

    boundary_mode = 'fixed'
    differentiation_method = '6th-Order'
    j = 0
    E = None
    V = 0
    x = None
    concentrations = None

    def _update(self):
        """Calculate the flux of chemical species."""
        self.set_E()
        self.dcdt = self.diffusive_flux() + \
            self.electromigration_flux() +\
            self.advection_flux()

    def set_current(self):
        """Calculate the current based on a fixed voltage drop."""
        self.state.current_density = (self.state.voltage /
                                      sum(self.dz / self.conductivity()))

    def set_E(self):
        """Calculate the electric field at each node."""
        self.state.field = -self.state.current_density/self.conductivity()

    def diffusive_flux(self):
        """Calculate flux due to diffusion."""
        cD = self.state.diffusivity * self.state.concentrations
        diffusion = \
            self.second_derivative(cD)
        return diffusion

    def electromigration_flux(self):
        """Calculate flux due to electromigration."""
        uc = self.state.mobility * self.state.concentrations
        electromigration = \
            self.first_derivative(uc * self.state.field)
        return electromigration

    def advection_flux(self):
        advection = (-self.state.bulk_flow *
                     self.first_derivative(self.state.concentrations))
        return advection

    def conductivity(self):
        """Calculate the conductivty at each location."""
        conductivity = np.sum(self.state.molar_conductivity
                              * self.state.concentrations, 0)
        return conductivity

    def node_flux(self):
        return np.zeros(self.state.nodes.shape)

    def pack(self, frame=None):
        if frame:
            return frame.concentrations.ravel()
        else:
            return self.dcdt.ravel()

    def unpack(self, packed):
        self.state.concentrations = \
            packed.reshape(self.state.concentrations.shape)
