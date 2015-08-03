def preconditioner(self):
    """Precondition the system by spacing the grid points."""
    # set up the interpolator to get the new parameters
    concentration_interpolator = interp1d(self.x,
                                          self.concentrations,
                                          kind='cubic')
    if self.area_variation:
        area_interpolator = interp1d(self.x,
                                     self.area,
                                     kind='cubic')

    # Update the flux calculator
    self.fluxer.update(self.x, self.area, self.concentrations)

    # Get the node cost from the flux calculator
    cost = self.fluxer.node_cost()
    cost = self.fluxer.differ.smooth(cost)

    # The last cost is poorly calculated, set it to an intermediate value
    cost[-1] = np.median(cost)

    # get the new grid parameters
    self.x = self._precondition_x(cost)
    self.concentrations = concentration_interpolator(self.x)
    if self.area_variation:
        self.area = area_interpolator(self.x)

    # equilibrate the new system.
    self.equilibrator.cH = self.equilibrator.pH = None
    self.equilibrator.equilibrate(self.concentrations)
    self.fluxer.update_ion_parameters(self.equilibrator)


def _precondition_x(self, cost):
    """Precondition the grid based on a cost."""
    new_x = np.cumsum(1/cost)
    new_x -= new_x[0]
    new_x *= max(self.x)/new_x[-1]
    return new_x
