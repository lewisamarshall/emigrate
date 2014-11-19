def set_current(self):
    """Calculate the current based on a fixed voltage drop."""
    if self.adaptive_grid is True:
        self.j = self.V/sum(self.dz * self.first_derivative(self.x) /
                            self.conductivity())
    else:
        self.j = self.V/sum(self.dz / self.conductivity())

def set_E(self):
    """Calculate the electric field at each node."""
    self.set_current()
    self.E = -self.j/self.conductivity()
    if self.adaptive_grid is True:
        self.E = self.E * self.first_derivative(self.x)

def set_derivatives(self):
    self.xz = self.first_derivative(self.x)
    self.xzz = self.second_derivative(self.x)

def reshaped_flux(self, t, state):
    """1-D flux function for ode solver."""
    self.t = t
    (self.x, self.concentrations) = self.decompose_state(state)
    self.set_derivatives()

    x_flux = self.node_flux()
    ion_flux = self.flux()
    flux = self.compose_state(x_flux, ion_flux)
    return flux

def set_current(self):
    """Calculate the current based on a fixed voltage drop."""
    if self.adaptive_grid is True:
        self.j = self.V/sum(self.dz * self.first_derivative(self.x) /
                            self.conductivity())
    else:
        self.j = self.V/sum(self.dz / self.conductivity())

def set_E(self):
    """Calculate the electric field at each node."""
    self.set_current()
    self.E = -self.j/self.conductivity()
    if self.adaptive_grid is True:
        self.E = self.E * self.first_derivative(self.x)

def flux(self):
    """Calculate the flux of chemical species."""
    self.set_E()
    total_flux = self.diffusive_flux() + \
        self.electromigration_flux() + \
        self.node_movement_flux()# + \
        # self.numerical_diffusion()
    total_flux = self.set_boundary(total_flux)
    return total_flux

def diffusive_flux(self):
    """Calculate flux due to diffusion."""
    cD = self.diffusivity * self.concentrations
    if self.adaptive_grid is True:
        diffusion = \
            (self.second_derivative(cD) -
             self.first_derivative(cD) * self.xzz / self.xz) / \
            (self.xz**2.)
    else:
        diffusion = \
            self.second_derivative(cD)
    return diffusion

def electromigration_flux(self):
    """Calculate flux due to electromigration."""
    uc = self.mobility * self.concentrations

    temp = (uc * (self.first_derivative(self.E) -(self.xzz/self.xz) * self.E)
            + self.first_derivative(uc) * self.E)/self.xz**2.
    if self.adaptive_grid is True:
        electromigration = temp
    else:
        electromigration = \
            self.first_derivative(uc * self.E)
        # print np.max(electromigration-temp)/np.max(electromigration)
    return electromigration

def numerical_diffusion(self):
    self.set_alpha()
    diff = self.alpha * (self.second_derivative(self.concentrations)-
                        self.first_derivative(self.first_derivative(self.concentrations)*
                         self.limit(self.first_derivative(self.concentrations))))
    return diff

def node_movement_flux(self):
    if self.adaptive_grid is True:
        node_movement = ((self.node_flux()-self.u) / self.xz) * \
            self.first_derivative(self.concentrations)
    else:
        node_movement = self.first_derivative(-self.u * self.concentrations)
    return node_movement

def set_boundary(self, flux):
    if self.boundary_mode == 'fixed':
        flux[:, 0] *= 0
        flux[:, -1] *= 0
    elif self.boundary_mode == 'characteristic':
        pass
    return flux

def node_flux(self):
    """Calculate the flux of nodes."""
    if self.adaptive_grid is True:
        flux = self.pointwave *\
            self.first_derivative(self.node_cost() *
                                  self.first_derivative(self.x))
        if False:
            window = np.bartlett(self.N_window)
            flux = np.convolve(flux, window, 'same')
        else:
            flux = self.differ.smooth(flux)
        flux[0, ] = flux[-1, ] = 0.
    else:
        flux = np.zeros(self.x.shape)
    return flux

def node_cost(self):
    """Calculate the cost function of each node."""
    deriv = np.abs(self.first_derivative(self.concentrations))
    cost = deriv / np.tile(np.nanmax(deriv, 1), (len(self.z), 1)).T
    cost = np.nanmax(cost, 0) + self.Kag
    return cost

def reshaped_flux(self, t, state):
    """1-D flux function for ode solver."""
    self.t = t
    (self.x, self.concentrations) = self.decompose_state(state)
    self.set_derivatives()

    x_flux = self.node_flux()
    ion_flux = self.flux()
    flux = self.compose_state(x_flux, ion_flux)
    return flux

def set_Kag(self):
    """Set the Kag parameter for spacing of low-gradient grid points."""
    self.Kag = ((self.N-self.NI)/self.NI) * self.Vthermal / self.V

def conductivity(self):
    """Calculate the conductivty at each location."""
    conductivity = np.sum(np.tile(self.molar_conductivity,
                                  (1, self.N))
                          * self.concentrations, 0)
    return conductivity
