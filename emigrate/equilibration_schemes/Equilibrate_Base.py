"""Define the class that other equilibration classes are derived from."""


class Equilibrate_Base(object):

    """A base class for ion equilibration of ions on grid points."""

    ions = None
    pH = 7
    concentrations = None
    mobility = None
    diffusivity = None
    molar_conductivity = None
    ionic_strength = None

    def __init__(self, ions, pH, concentrations):
        """Initialize with ions, pH concentrations."""
        self.ions = ions
        self.pH = pH
        self.concentrations = concentrations
        self.set_arrays()

    def calc_equilibrium(self):
        """Call subroutines for mobility, diffusivity, and conductivity."""
        self.mobility = self.concentrations * 0.
        self.diffusivity = self.concentrations * 0.
        self.molar_conductivity = self.concentrations * 0.

    def equilibrate(self, concentrations):
        """Perform full equlibration routine."""
        self.concentrations = concentrations
        self.calc_equilibrium()
        return (self.mobility, self.diffusivity, self.molar_conductivity)

    def set_arrays(self):
        """Placeholder for initialization setup."""
        pass
