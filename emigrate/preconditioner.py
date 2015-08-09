from scipy.interpolate import interp1d
import numpy as np


def preconditioner(state, fluxer):
    """Precondition the system by spacing the grid points."""
    # set up the interpolator to get the new parameters
    concentration_interpolator = interp1d(state.nodes,
                                          state.concentrations,
                                          kind='cubic')

    # Get the node cost from the flux calculator
    fluxer.update()
    cost = fluxer.node_cost()
    cost = fluxer.differ.smooth(cost)

    # The last cost is poorly calculated, set it to an intermediate value
    cost[-1] = np.median(cost)

    # get the new grid parameters
    state.nodes = _precondition_x(state, cost)
    state.concentrations = concentration_interpolator(state.nodes)
    if np.size(state.area) > 1:
        area_interpolator = interp1d(state.nodes,
                                     state.area,
                                     kind='cubic')
        state.area = area_interpolator(state.nodes)


def _precondition_x(state, cost):
    """Precondition the grid based on a cost."""
    new_x = np.cumsum(1/cost)
    new_x -= new_x[0]
    new_x *= max(state.nodes)/new_x[-1]
    return new_x
