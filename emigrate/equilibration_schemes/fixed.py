import numpy as np
from Equilibrate_Base import Equilibrate_Base

class Fixed(Equilibrate_Base):

    def calc_equilibrium(self):
        self.calc_diffusivity()
        self.calc_mobility()
        self.calc_molar_conductivity()

    def calc_diffusivity(self):
        self.diffusivity = np.array([[ion.diffusivity(self.pH)]
                                    for ion in self.ions])

    def calc_mobility(self):
        self.mobility = np.array([[ion.effective_mobility(self.pH)]
                                  for ion in self.ions])

    def calc_molar_conductivity(self):
        self.molar_conductivity = np.array([[ion.molar_conductivity(self.pH)]
                                            for ion in self.ions])
