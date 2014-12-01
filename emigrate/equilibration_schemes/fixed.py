"""An equilibration class for fixed pH problems."""
import numpy as np
from Equilibrate_Base import Equilibrate_Base
# pylint: disable = W0232, E1101, W0201


class Fixed(Equilibrate_Base):

    """An equilibration class for fixed pH problems."""

    def calc_equilibrium(self):
        """Calculate the equilibrium properties."""
        self.calc_diffusivity()
        self.calc_mobility()
        self.calc_molar_conductivity()

    def calc_diffusivity(self):
        """Calculate the diffusivity."""
        self.diffusivity = np.array([[ion.diffusivity(self.pH)]
                                    for ion in self.ions])

    def calc_mobility(self):
        """Calculate the effective mobility."""
        self.mobility = np.array([[ion.effective_mobility(self.pH)]
                                  for ion in self.ions])

    def calc_molar_conductivity(self):
        """Calculate the conductivity."""
        self.molar_conductivity = np.array([[ion.molar_conductivity(self.pH)]
                                            for ion in self.ions])
