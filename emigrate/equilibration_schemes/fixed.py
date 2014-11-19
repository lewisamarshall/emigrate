import numpy as np

class Fixed(object):
    def __init__(self, ions, pH):
        self.ions = ions
        self.pH = pH

    def set_ion_properties(self):
        """Set the properties of ions in the system."""
        self.diffusivity = np.array([[ion.diffusivity(self.pH)]
                                    for ion in self.ions])
        self.mobility = np.array([[ion.effective_mobility(self.pH)]
                                  for ion in self.ions])
        self.molar_conductivity = np.array([[ion.molar_conductivity(self.pH)]
                                            for ion in self.ions])
        return (self.mobility, self.diffusivity, self.molar_conductivity)
