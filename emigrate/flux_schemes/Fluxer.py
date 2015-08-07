from Differentiate import Differentiate
import numpy as np
import warnings

# #TODO:40 Fix boundary characteristics
# from boundary_characteristic import (boundary_characteristic,
#                                      _get_characteristic_matricies,
#                                      _a_matrix)


class Fluxer(object):

    """A base class for the flux modules used by emigrate."""

    # Differentiation information
    differentiation_method = '6th-Order'

    # Boundary condition information
    boundary_mode = 'characteristic'
    nonnegative = True

    # Field info
    mode = 'voltage'

    # Dependant state information
    # #TODO:0 Change implimentation
    _area = 1.

    # Reference Frame
    frame = None
    edge = 'right'
    _frame_velocity = 0

    def __init__(self, state):
        """Initialize the compact flux solver."""
        self.state = state

        # Prepare the grid points from the state nodes
        self._prep_domain()

        # Prepare the voltage/current and bulk flow mode from the state.
        self.current_density = self.state.current_density
        if self.state.voltage:
            self.mode = 'voltage'
            if self.current_density:
                warnings.warn(
                    'state has both current and voltage. Using voltage.'
                    )
        else:
            self.mode = 'current'

        # use state area if it exists, otherwise default to _area
        if self.state.area is not None:
            self._area = self.state.area
        self._area = np.array(self._area)
        if self._area.size > 1:
            self.area_variation = True

        # Create the differentiation state.
        self.differ = Differentiate(self.N, self.dz,
                                    method=self.differentiation_method)

    def _prep_domain(self):
        self.z = np.linspace(min(self.state.nodes),
                             max(self.state.nodes),
                             len(self.state.nodes))
        self.N = self.state.nodes.size
        self.dz = self.z[1]-self.z[0]

    def first_derivative(self, x_input):
        """Calculate the first derivative with respect to z."""
        return self.differ.first_derivative(x_input.T).T

    def second_derivative(self, x_input):
        """Calculate the second derivative with respect to z."""
        return self.differ.second_derivative(x_input.T).T

    def set_boundary(self):
        """Set the boundary condition at the domain edges."""
        if self.boundary_mode == 'fixed':
            self.dcdt[:, 0] = self.dcdt[:, -1] = 0.
        elif self.boundary_mode == 'characteristic':
            self.dcdt[:, 0] = self.boundary_characteristic('left')
            self.dcdt[:, -1] = self.boundary_characteristic('right')

    def update(self):
        self._update()
        self.set_boundary()

        # Impose nonnegativity constraint.
        if self.nonnegative is True:
            self.dcdt = self._impose_nonnegativity(self.state.concentrations,
                                                   self.dcdt)

        # Update the reference frame for the next time step.
        if self.frame is not None:
            self._update_reference_frame()

    def _impose_nonnegativity(self, concentrations, dcdt):
        dcdt = np.where(np.greater(concentrations, 0),
                        dcdt,
                        np.maximum(0, dcdt))
        return dcdt

    def _update_reference_frame(self):
        if self.edge == 'right':
            E = self.state.field[-1]
            pH = self.state.pH[-1]
        elif self.edge == 'left':
            E = self.state.field[0]
            pH = self.state.pH[0]
        else:
            raise RuntimeError('Edge must be left or right.')
        self._frame_velocity = E * self.frame.effective_mobility(pH)

    def _update(self):
        raise NotImplementedError

    def pack(self, frame):
        raise NotImplementedError

    def unpack(self, packed, frame):
        raise NotImplementedError
