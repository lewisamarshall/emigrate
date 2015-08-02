"""Define the class that other Equilibrator classes are derived from."""


class Equilibrator(object):

    """A base class for ion equilibration of ions on grid points."""

    # Initialization variables.
    ions = None
    concentrations = None

    # Required core properties
    # TODO: Remove and place in state
    pH = None
    cH = None
    ionic_strength = None

    # Required secondary properties
    # These are calculated from the core properties after
    # equilibration is achieved.
    # TODO: remove and place in state.
    mobility = None
    diffusivity = None
    molar_conductivity = None

    def __init__(self, state):
        """Initialize with ions, pH concentrations."""
        self.state = state

    def equilibrate(self):
        """Equilibrate a new concentration set.

        This is the only public method in an Equilibrator.
        """
        raise NotImplementedError
