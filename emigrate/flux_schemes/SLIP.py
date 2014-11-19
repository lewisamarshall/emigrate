def set_alpha(self):
    """Set alpha for dissipation."""
    self.set_characteristic()
    self.alpha = 0.5 * np.maximum(np.fabs(self.characteristic/self.tile_m(self.xz)),0)

def set_characteristic(self):
    """Calculate the characteristic speed of paramters."""
    self.characteristic = self.u + self.tile_m(self.E)*self.tile_n(self.mobility) -\
        self.node_flux()

def limit(self, x_input):
    return self.limiter.limit(x_input)
