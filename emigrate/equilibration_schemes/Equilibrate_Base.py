class Equilibrate_Base(object):
    """A base class for ion equilibration."""
    ions = None
    pH = 7
    concentrations = None
    mobility = None
    diffusivity = None
    molar_conductivity = None

    def __init__(self, ions, pH, concentrations):
        self.ions = ions
        self.pH = pH
        self.concentrations = concentrations

    def calc_equilibrium(self):
        self.mobility = self.concentrations * 0.
        self.diffusivity = self.concentrations * 0.
        self.molar_conductivity = self.concentrations * 0.

    def equilibrate(self, concentrations):
        self.concentrations = concentrations
        self.calc_equilibrium()
        return (self.mobility, self.diffusivity, self.molar_conductivity)
