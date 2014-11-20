from Differentiate import Differentiate


class _Flux_Base(object):

    """A base class for the flux modules used by emigrate."""

    use_adaptive_grid = False
    boundary_mode = 'characteristic'
    differentiation_method = '6th-Order'
    j = 0
    E = None
    V = 0
    x = None
    concentrations = None

    def __init__(self, N, dz, V,  mobility, diffusivity, molar_conductivity):
        """Initialize the compact flux solver."""
        self.N = N
        self.dz = dz
        self.V = V
        self.mobility = mobility
        self.diffusivity = diffusivity
        self.molar_conductivity = molar_conductivity
        self.differ = Differentiate(N, dz, method=self.differentiation_method)

    def first_derivative(self, x_input):
        """Calculate the first derivative with respect to z."""
        return self.differ.first_derivative(x_input.T).T

    def second_derivative(self, x_input):
        """Calculate the second derivative with respect to z."""
        return self.differ.second_derivative(x_input.T).T

    def set_boundary(self, flux):
        """Set the boundary condition at the domain edges."""
        if self.boundary_mode == 'fixed':
            flux[:, 0] *= 0
            flux[:, -1] *= 0
        elif self.boundary_mode == 'characteristic':
            pass
        return flux

    def flux(self, x, concentrations):
        """Calculate the flux of chemical species."""
        self.x = x
        self.concentrations = concentrations
        return self.concentrations * 0.
