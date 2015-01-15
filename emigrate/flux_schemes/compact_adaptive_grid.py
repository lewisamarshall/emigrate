from compact import Compact
import numpy as np

class CompactAdaptive(Compact):
    pointwave = 1

    def flux(self, x, concentrations):
        """Caclulate the flux of chemical species."""
        self.x = x
        self.concentrations = concentrations
        self.set_derivatives()
        self.set_E()

        total_flux = (self.diffusive_flux() +
                      self.electromigration_flux() +
                      self.node_movement_flux())

        total_flux = self.set_boundry(total_flux)
        return total_flux

    def set_derivatives(self):
        self.xz = self.first_derivative(self.x)
        self.xzz = self.second_derivative(self.x)

    def set_E(self):
        """Calculate the electric field at each node."""
        self.set_current()
        self.E = -self.j/self.conductivity()
        # if self.adaptive_grid is True:
        self.E = self.E * self.first_derivative(self.x)

    def set_current(self):
        """Calculate the current based on a fixed voltage drop."""
        self.j = self.V/sum(self.dz * self.first_derivative(self.x) /
                            self.conductivity())

    # def conductivity(self):
    #     """Calculate the conductivty at each location."""
    #     conductivity = np.sum(np.tile(self.molar_conductivity,
    #                                   (1, self.N))
    #                           * self.concentrations, 0)

    def diffusive_flux(self):
        """Calculate flux due to diffusion."""
        cD = self.diffusivity * self.concentrations
        diffusion = \
            (self.second_derivative(cD) -
             self.first_derivative(cD) * self.xzz / self.xz) / \
            (self.xz**2.)
        return diffusion

    def electromigration_flux(self):
        """Calculate flux due to electromigration."""
        uc = self.mobility * self.concentrations

        electromigration = (uc * (self.first_derivative(self.E) -
                            (self.xzz/self.xz) * self.E)
                            + self.first_derivative(uc) * self.E) / self.xz**2.
        return electromigration

    def node_movement_flux(self):
        node_movement = ((self.node_flux()-self.u) / self.xz) * \
            self.first_derivative(self.concentrations)
        return node_movement

    def set_boundary(self, flux):
        if self.boundary_mode == 'fixed':
            flux[:, 0] *= 0
            flux[:, -1] *= 0
        elif self.boundary_mode == 'characteristic':
            pass
        return flux

    def node_flux(self):
        """Calculate the flux of nodes."""
        flux = self.pointwave *\
            self.first_derivative(self.node_cost() *
                                  self.first_derivative(self.x))
        if False:
            window = np.bartlett(self.N_window)
            flux = np.convolve(flux, window, 'same')
        else:
            flux = self.differ.smooth(flux)
        flux[0, ] = flux[-1, ] = 0.

        return flux

    def node_cost(self):
        """Calculate the cost function of each node."""
        deriv = np.abs(self.first_derivative(self.concentrations))
        cost = deriv / np.tile(np.nanmax(deriv, 1), (len(self.z), 1)).T
        cost = np.nanmax(cost, 0) + self.Kag
        return cost

    def set_Kag(self):
        """Set the Kag parameter for spacing of low-gradient grid points."""
        self.Kag = ((self.N-self.NI)/self.NI) * self.Vthermal / self.V

        return conductivity
