"""An equilibration class for fixed pH problems."""
import numpy as np
from Equilibrator import Equilibrator
# pylint: disable = W0232, E1101, W0201


class Fixed(Equilibrator):

    """An equilibration class for fixed pH problems."""

    def __init__(self, state):
        super(Fixed, self).__init__(state)
        if self.state.pH is None:
            self.state.pH = 7.
        if self.state.cH is None:
            self.state.cH = 10**(-self.state.pH)
        self.state.water_diffusive_conductivity = 0
        self.state.water_conductivity = 0

    def equilibrate(self):
        """Calculate the equilibrium properties."""
        self._calc_diffusivity()
        self._calc_mobility()
        self._calc_molar_conductivity()

    def _calc_diffusivity(self):
        """Calculate the diffusivity."""
        self.state.diffusivity = np.array([[ion.diffusivity(self.state.pH)]
                                          for ion in self.state.ions])

    def _calc_mobility(self):
        """Calculate the effective mobility."""
        self.state.mobility = np.array([[ion.effective_mobility(self.state.pH)]
                                       for ion in self.state.ions])

    def _calc_molar_conductivity(self):
        """Calculate the conductivity."""
        self.state.molar_conductivity = \
            np.array([[ion.molar_conductivity(self.state.pH)]
                     for ion in self.state.ions])
