"""Defines an equlibration scheme with pH calculation."""

import numpy as np
from Equilibrate_Base import Equilibrate_Base
from math import log10
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
    z0_matrix = None
    ionization_fraction = None
    absolute_mobility = None
    _Lpm3 = 1000
    _F = 96485.34
    _kB = 8.617e-6        # EV/K
    T = 25

    def set_arrays(self):
        """Prepare arrays to solve problems during initialization."""
        self.set_l_matrix()
        self.set_Q()
        self.set_Pmat()
        self.set_z0_matrix()
        self.set_absolute_mobility()

    def set_absolute_mobility(self):
        """Build the absolute mobility matrix."""
        self.absolute_mobility = np.array([i.absolute_mobility
                                           for i in self.ions])

    def set_l_matrix(self):
        """Build the L matrix."""
        # Find the order of the polynomial. This is the maximum
        # size of the list of charge states in an ion.
        self.max_columns = max([max(i.z)-min(i.z)+2 for i in self.ions])
        self.m_ions = len(self.ions)
        self.nodes = self.concentrations.shape[1]

        # Set up the matrix of Ls, the multiplication
        # of acidity coefficients for each ion.
        self.l_matrix = np.array([np.resize(i.L(I=0),
                                 [self.max_columns]) for i in self.ions])

    def set_z0_matrix(self):
        """Build the matrix of valences."""
        self.z0_matrix = np.array([i.z0 for i in self.ions])

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
            z_list = np.resize(self.ions[i].z0, [self.max_columns])

            Mmod = self.l_matrix.copy()
            Mmod[i, :] *= np.array(z_list)

            Pi = 1.
            for k in range(self.m_ions):
                Pi = np.convolve(Pi, Mmod[k, :])

            Pi = np.convolve([0.0, 1.0], Pi)  # Convolve with P2
            self.PMat.append(Pi)
        self.PMat = np.array(self.PMat, ndmin=2)[:, :, np.newaxis]

    def calc_equilibrium(self):
        """Calculate equilibrium."""
        self.calc_pH()
        self.calc_ionization_fraction()
        self.calc_mobility()
        self.calc_diffusivity()
        self.calc_molar_conductivity()

    def calc_pH(self):
        """Return the pH of the object."""
        # Construct P matrix

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
        # reverse order for poly function
        self.pH = []
        self.cH = []
        for i in range(self.nodes):
            cH = np.roots(poly[:, i])

            cH = [c for c in cH if c.real > 0 and c.imag == 0]
            if cH:
                cH = float(cH[-1].real)
            else:
                print 'Failed to find pH.'

        # Convert to pH. Use the activity to correct the calculation.
            self.pH.append(-log10(cH))
            self.cH.append(cH)

    def calc_mobility(self):
        """Calculate effective mobility."""
        self.mobility = np.sum(self.ionization_fraction *
                               self.absolute_mobility[:, :, np.newaxis], 1)

    def calc_diffusivity(self):
        """Calculate diffusivity."""
        self.diffusivity = (self.absolute_mobility[:, :, np.newaxis] *
                            self.ionization_fraction /
                            (self.z0_matrix[:, :, np.newaxis]+.5)) *\
            self._kB * (self.T + 273.15)
        self.diffusivity = np.sum(self.diffusivity, 1)

    def calc_molar_conductivity(self):
        """Calculate molar conductivity."""
        self.molar_conductivity = self._Lpm3 * self._F * \
            np.sum(self.z0_matrix[:, :, np.newaxis] *
                   self.ionization_fraction *
                   self.absolute_mobility[:, :, np.newaxis], 1)

    def calc_ionization_fraction(self):
        """Calculate ionization fraction."""
        # Get the vector of products of acidity constants.
        # Compute the concentration of H+ from the pH.
        # Calculate the numerator of the function for ionization fraction.
        i_frac_vector = self.l_matrix[:, :, np.newaxis] *\
            self.cH**self.z0_matrix[:, :, np.newaxis]

        # Calculate the vector of ionization fractions
        # Filter out the neutral fraction
        denom = np.sum(i_frac_vector, 1)

        self.ionization_fraction = i_frac_vector/denom[:, np.newaxis, :]
