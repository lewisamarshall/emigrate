"""Define the class that other Equilibrator classes are derived from."""


class Equilibrator(object):

    """A base class for ion equilibration of ions on grid points."""

    # Initialization variables.
    ions = None
    concentrations = None

    # Required core properties
    pH = None
    cH = None
    ionic_strength = None

    # Required secondary properties
    # These are calculated from the core properties after
    # equilibration is achieved.
    mobility = None
    diffusivity = None
    molar_conductivity = None

    def __init__(self, ions, concentrations):
        """Initialize with ions, pH concentrations."""
        self.ions = ions
        self.concentrations = concentrations

    def equilibrate(self, concentrations):
        """Equilibrate a new concentration set.

        This is the only public method in an Equilibrator.
        """
        self.concentrations = concentrations
        self._equilibrate()

    def _equilibrate(self):
        """Placeholder equilibration scheme."""
        raise NotImplementedError
