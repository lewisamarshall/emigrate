"""Defines an equlibration scheme with pH calculation."""

import numpy as np
from Equilibrator import Equilibrator
from Multiroot import Multiroot
# pylint: disable=W0232, E1101, W0201, E1103

# Physical Constants
k_water = 1E-14
lpm3 = 1000
faraday = 96485.34          # Faraday's const.[C/mol]
boltzmann = 8.617e-6        # EV/K
temperature = 25
temperature_K = temperature + 273.15
h_mobility = 362E-9/faraday
oh_mobility = 205E-9/faraday
h_diffusivity = h_mobility / 1 * boltzmann * (temperature_K)
oh_diffusivity = oh_mobility / -1 * boltzmann * (temperature_K)


class VariablepH(Equilibrator):

    """A solver that uses calculates pH at each time point."""

    # Private properties used during calculations
    _l_matrix = None
    _Q = None
    _PMat = None
    _z0 = None
    _z = None
    _index_0 = None

    # New public properties
    # TODO: move to state.
    ionization_fraction = None
    absolute_mobility = None
    water_conductivity = None
    water_diffusive_conductivity = None

    def __init__(self, state):
        super(VariablepH, self).__init__(state)
        self._prepare_arrays()
        self._multiroot = Multiroot()

    def _prepare_arrays(self):
        """Prepare arrays to solve problems during initialization."""
        self._set_z_index()
        self._set_l_matrix()
        self._set_Q()
        self._set_Pmat()
        self._set_absolute_mobility()

    def equilibrate(self):
        """Calculate equilibrium."""
        self._calc_pH()
        self._calc_ionization_fraction()
        self._calc_mobility()
        self._calc_diffusivity()
        self._calc_molar_conductivity()
        self._calc_water_conductivity()
        self._calc_water_diffusive_conductivity()

    def _set_z_index(self):
        """Set the valence indices."""
        all_z = []
        for i in self.state.ions:
            all_z.extend(i.z0)
        self._z0 = range(min(all_z), max(all_z)+1)
        self._index_0 = self._z0.index(0)

        self._z0 = np.array(self._z0)
        self._z = self._z0.tolist()
        self._z.pop(self._index_0)
        self._z = np.array(self._z)

    def _align_zero(self, value, z0):
        """Align ion properties with the zero of the matrix."""
        local_index = z0.index(0)
        local_len = len(z0)
        pre_pad = self._index_0 - local_index
        post_pad = len(self._z0) - local_len - pre_pad
        return np.pad(value,
                      (pre_pad, post_pad),
                      'constant', constant_values=(0))

    def _set_absolute_mobility(self):
        """Build the absolute mobility matrix."""
        self.absolute_mobility = []
        for i in self.state.ions:
            self.absolute_mobility.append(self._align_zero(i.absolute_mobility,
                                                           i.z0))
        self.absolute_mobility = np.array(self.absolute_mobility)

    def _set_l_matrix(self):
        """Build the L matrix."""
        # Set up the matrix of Ls, the multiplication
        # of acidity coefficients for each ion.
        self._l_matrix = []
        for i in self.state.ions:
            self._l_matrix.append(self._align_zero(i.L(I=0), i.z0))
        self._l_matrix = np.array(self._l_matrix)

    def _set_Q(self):
        """Build the Q matrix for pH solving."""
        # Construct Q vector.
        self._Q = 1.
        for j in range(len(self.state.ions)):
            self._Q = np.convolve(self._Q, self._l_matrix[j, :])

        # Convolve with water dissociation.
        self._Q = np.convolve(self._Q, [-k_water, 0.0, 1.0])

    def _set_Pmat(self):
        """Build the Pmat Matrix for pH solving."""
        self._PMat = []
        for i in range(len(self.state.ions)):

            Mmod = self._l_matrix.copy()
            Mmod[i, :] *= self._z0

            Pi = 1.
            for k in range(len(self.state.ions)):
                Pi = np.convolve(Pi, Mmod[k, :])

            Pi = np.convolve([0.0, 1.0], Pi)  # Convolve with P2
            self._PMat.append(Pi)
        self._PMat = np.array(self._PMat, ndmin=2)[:, :, np.newaxis]

    def _calc_pH(self):
        """Return the pH of the object."""
        # Multiply P matrix by concentrations, and sum.
        P = np.sum(self._PMat *
                   np.array(self.state.concentrations)[:, np.newaxis, :], 0)

        # Construct polynomial. Change the shapes, then reverse  order
        if P.shape[0] < self._Q.shape[0]:
            P.resize((self._Q.shape[0], P.shape[1]))
        elif P.shape[0] > self._Q.shape[0]:
            self._Q.resize(P.shape[0])
        poly = (P + self._Q[:, np.newaxis])[::-1]

        self.cH = self.state.cH = self._multiroot(poly, self.cH)

        self.pH = self.state.pH = -np.log10(self.cH)

        if any(np.isnan(self.pH)):
            print 'pH:', self.pH
            print 'cH:', self.cH
            raise RuntimeError("Couldn't find correct pH.")

    def _calc_mobility(self):
        """Calculate effective mobility."""
        self.mobility = np.sum(self.ionization_fraction *
                               self.absolute_mobility[:, :, np.newaxis], 1)

    def _calc_diffusivity(self):
        """Calculate diffusivity."""
        self.diffusivity = (self.absolute_mobility[:, :, np.newaxis] *
                            self.ionization_fraction /
                            (self._z[np.newaxis, :, np.newaxis])) *\
            boltzmann * (temperature_K)
        self.diffusivity = np.sum(self.diffusivity, 1)

    def _calc_molar_conductivity(self):
        """Calculate molar conductivity."""
        self.molar_conductivity = lpm3 * faraday * \
            np.sum(self._z[np.newaxis, :, np.newaxis] *
                   self.ionization_fraction *
                   self.absolute_mobility[:, :, np.newaxis], 1)

    def _calc_ionization_fraction(self):
        """Calculate ionization fraction."""
        # Calculate the numerator of the function for ionization fraction.
        i_frac_vector = self._l_matrix[:, :, np.newaxis] *\
            self.cH**self._z0[np.newaxis, :, np.newaxis]

        # Calculate the vector of ionization fractions
        denom = np.sum(i_frac_vector, 1)

        # Filter out the uncharged state.
        self.ionization_fraction = i_frac_vector/denom[:, np.newaxis, :]
        self.ionization_fraction = np.delete(self.ionization_fraction,
                                             self._index_0,
                                             axis=1)

    def _calc_water_conductivity(self):
        self.water_conductivity = (self.cH * h_mobility +
                                   k_water / self.cH * oh_mobility)

    def _calc_water_diffusive_conductivity(self):
        self.water_diffusive_conductivity = \
            (self.cH * h_diffusivity -
             k_water/self.cH * oh_diffusivity) * faraday
