"""Define the class that other Equilibrator classes are derived from."""


class Equilibrator(object):

    """A base class for ion equilibration of ions on grid points."""

    def __init__(self, state):
        """Initialize with ions, pH concentrations."""
        self.state = state

    def equilibrate(self):
        """Equilibrate a new concentration set.

        This is the only public method in an Equilibrator.
        """
        raise NotImplementedError
