"""Defines an equlibration scheme with pH calculation."""

import numpy as np
from Equilibrate_Base import Equilibrate_Base
from multiroot import multiroot
# pylint: disable=W0232, E1101, W0201, E1103


class Variable_pH(Equilibrate_Base):

    """A solver that uses calculates pH at each time point."""

    Kw = 1E-14
    l_matrix = None
    Q = None
    max_columns = None
    m_ions = None
    nodes = None
    PMat = None
    pH = None
    cH = None
    ionization_fraction = None
    absolute_mobility = None
    _Lpm3 = 1000
    _F = 96485.34
    _kB = 8.617e-6        # EV/K
    T = 25
    z0 = None
    index_0 = None
    F = 9.65E4       # Faraday's const.[C/mol]
    h_mobility = 362E-9/F
    oh_mobility = 205E-9/F

    h_diffusivity = h_mobility / 1 * _kB * (T+273.15)
    oh_diffusivity = oh_mobility / -1 * _kB * (T+273.15)
    water_conductivity = None

    def set_arrays(self):
        """Prepare arrays to solve problems during initialization."""
        self.set_z_index()
        self.set_l_matrix()
        self.set_Q()
        self.set_Pmat()
        self.set_absolute_mobility()

    def calc_equilibrium(self):
        """Calculate equilibrium."""
        self.calc_pH()
        self.calc_ionization_fraction()
        self.calc_mobility()
        self.calc_diffusivity()
        self.calc_molar_conductivity()
        self.calc_water_conductivity()
        self.calc_water_diffusive_conductivity()

    def set_z_index(self):
        """Set the valence indices."""
        all_z = []
        for i in self.ions:
            all_z.extend(i.z0)
        all_z = set(all_z)
        self.z0 = range(min(all_z), max(all_z)+1)
        self.index_0 = self.z0.index(0)

        self.m_ions = len(self.ions)
        self.nodes = self.concentrations.shape[1]
        self.max_columns = len(self.z0)
        self.z0 = np.array(self.z0)
        self.z = self.z0.tolist()
        self.z.pop(self.index_0)
        self.z = np.array(self.z)

    def align_zero(self, value, z0):
        """Align ion properties with the zero of the matrix."""
        local_index = z0.index(0)
        local_len = len(z0)
        pre_pad = self.index_0 - local_index
        post_pad = len(self.z0) - local_len - pre_pad
        return np.pad(value,
                      (pre_pad, post_pad),
                      'constant', constant_values=(0))

    def set_absolute_mobility(self):
        """Build the absolute mobility matrix."""
        self.absolute_mobility = []
        for i in self.ions:
            self.absolute_mobility.append(self.align_zero(i.absolute_mobility,
                                                          i.z0))
        self.absolute_mobility = np.array(self.absolute_mobility)

    def set_l_matrix(self):
        """Build the L matrix."""
        # Set up the matrix of Ls, the multiplication
        # of acidity coefficients for each ion.
        self.l_matrix = []
        for i in self.ions:
            self.l_matrix.append(self.align_zero(i.L(I=0), i.z0))
        self.l_matrix = np.array(self.l_matrix)

    def set_Q(self):
        """Build the Q matrix for pH solving."""
        # Construct Q vector.
        self.Q = 1.
        for j in range(self.m_ions):
            self.Q = np.convolve(self.Q, self.l_matrix[j, :])

        # Convolve with water dissociation.
        self.Q = np.convolve(self.Q, [-self.Kw, 0.0, 1.0])

    def set_Pmat(self):
        """Build the Pmat Matrix for pH solving."""
        self.PMat = []
        for i in range(self.m_ions):

            Mmod = self.l_matrix.copy()
            Mmod[i, :] *= self.z0

            Pi = 1.
            for k in range(self.m_ions):
                Pi = np.convolve(Pi, Mmod[k, :])

            Pi = np.convolve([0.0, 1.0], Pi)  # Convolve with P2
            self.PMat.append(Pi)
        self.PMat = np.array(self.PMat, ndmin=2)[:, :, np.newaxis]

    def calc_pH(self):
        """Return the pH of the object."""
        # Multiply P matrix by concentrations, and sum.
        P = np.sum(self.PMat *
                   np.array(self.concentrations)[:, np.newaxis, :], 0)

        # Construct polynomial. Change the shapes, then reverse  order
        if P.shape[0] < self.Q.shape[0]:
            P.resize((self.Q.shape[0], P.shape[1]))
        elif P.shape[0] > self.Q.shape[0]:
            self.Q.resize(P.shape[0])
        poly = (P+self.Q[:, np.newaxis])[::-1]

        # Solve Polynomial for concentration
        if self.pH is None or self.cH is None:
            self.cH = []

            for i in range(self.nodes):
                cH = np.roots(poly[:, i])

                cH = [c for c in cH if c.real > 0 and c.imag == 0]
                if cH:
                    cH = float(cH[-1].real)
                else:
                    print 'Failed to find pH.'

            # Convert to pH. Use the activity to correct the calculation.
                self.cH.append(cH)

            # convert to numpy array before ending
            self.cH = np.array(self.cH)

        else:
            self.cH = multiroot(poly, self.cH)


        self.pH = -np.log10(self.cH)

    def calc_mobility(self):
        """Calculate effective mobility."""
        self.mobility = np.sum(self.ionization_fraction *
                               self.absolute_mobility[:, :, np.newaxis], 1)

    def calc_diffusivity(self):
        """Calculate diffusivity."""
        self.diffusivity = (self.absolute_mobility[:, :, np.newaxis] *
                            self.ionization_fraction /
                            (self.z[np.newaxis, :, np.newaxis])) *\
            self._kB * (self.T + 273.15)
        self.diffusivity = np.sum(self.diffusivity, 1)

    def calc_molar_conductivity(self):
        """Calculate molar conductivity."""
        self.molar_conductivity = self._Lpm3 * self._F * \
            np.sum(self.z[np.newaxis, :, np.newaxis] *
                   self.ionization_fraction *
                   self.absolute_mobility[:, :, np.newaxis], 1)

    def calc_ionization_fraction(self):
        """Calculate ionization fraction."""
        # Calculate the numerator of the function for ionization fraction.
        i_frac_vector = self.l_matrix[:, :, np.newaxis] *\
            self.cH**self.z0[np.newaxis, :, np.newaxis]

        # Calculate the vector of ionization fractions
        denom = np.sum(i_frac_vector, 1)

        # Filter out the uncharged state.

        self.ionization_fraction = i_frac_vector/denom[:, np.newaxis, :]
        self.ionization_fraction = np.delete(self.ionization_fraction,
                                             self.index_0,
                                             axis=1)

    def calc_water_conductivity(self):
        self.water_conductivity = \
            self.cH * self.h_mobility + \
            self.Kw/self.cH * self.oh_mobility


    def calc_water_diffusive_conductivity(self):
        self.water_diffusive_conductivity = \
            (self.cH * self.h_diffusivity - \
            self.Kw/self.cH * self.oh_diffusivity)*self.F
