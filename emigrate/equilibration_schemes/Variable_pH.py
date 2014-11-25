import numpy as np
from Equilibrate_Base import Equilibrate_Base
import scipy
from math import log10
# pylint: disable=W0232, E1101, W0201


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

    def set_arrays(self):
        self.set_l_matrix()
        self.set_Q()
        self.set_Pmat()

    def set_l_matrix(self):
        # Find the order of the polynomial. This is the maximum
        # size of the list of charge states in an ion.
        self.max_columns = max([max(i.z)-min(i.z)+2 for i in self.ions])
        self.m_ions = len(self.ions)
        self.nodes = self.concentrations.shape[1]

        # Set up the matrix of Ls, the multiplication
        # of acidity coefficients for each ion.
        self.l_matrix = np.array([np.resize(i.L(I=0),
                             [self.max_columns]) for i in self.ions])

    def set_Q(self):
        # Construct Q vector.
        self.Q = 1.
        for j in range(self.m_ions):
            self.Q = np.convolve(self.Q, self.l_matrix[j, :])

        # Convolve with water dissociation.
        self.Q = np.convolve(self.Q, [-self.Kw, 0.0, 1.0])

    def set_Pmat(self):
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
        self.calc_mobility()
        self.calc_diffusivity()
        self.calc_molar_conductivity()

    def calc_pH(self):
        """Return the pH of the object."""
        # Construct P matrix

        # Multiply P matrix by concentrations, and sum.
        P = np.sum(self.PMat *
                   np.array(self.concentrations)[:, np.newaxis,:], 0)

        # Construct polynomial. Change the shapes as needed, then reverse  order
        if P.shape[0] < self.Q.shape[0]:
            P.resize((self.Q.shape[0], P.shape[1]))
        elif P.shape[0] > self.Q.shape[0]:
            Q.resize(P.shape[0])
        poly = (P+self.Q[:,np.newaxis])[::-1]
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
        self.mobility = np.array([[ion.effective_mobility(self.pH[j]) for j in range(self.concentrations.shape[1])]
                for ion in self.ions])

        # effective_mobility = sum([f*m for (f, m) in zip(i_frac, actual_mobility)])


    def calc_diffusivity(self):
        self.diffusivity =  np.array([[ion.diffusivity(self.pH[j]) for j in range(self.concentrations.shape[1])]
                for ion in self.ions])

    #     diffusivity = sum([m * f / float(z) for
    #                   m,f,z in zip(self.absolute_mobility,
    #                                self.ionization_fraction(pH),
    #                                self.z
    #                                )]) * self._kB * (self.T + 273.15)
    # diffusivity/=sum(self.ionization_fraction(pH))

    def calc_molar_conductivity(self):
        self.molar_conductivity =  np.array([[ion.molar_conductivity(self.pH[j]) for j in range(self.concentrations.shape[1])]
                for ion in self.ions])

        # m_conductivity = (self._Lpm3 * self._F *
        #               sum(z * f * m for (z, f, m)
        #                   in zip(self.z, i_frac, actual_mobility)))

    def calc_ionization_fraction(self):
        # Get the vector of products of acidity constants.
        L = self.L(I)
        # Compute the concentration of H+ from the pH.
        cH = 10**(-pH)/self.activity_coefficient(I, [1])[0]

        # Calculate the numerator of the function for ionization fraction.
        i_frac_vector = [Lp * cH ** z for (Lp, z) in zip(L, self.z0)]

        # Calculate the vector of ionization fractions
        # Filter out the neutral fraction
        denom = sum(i_frac_vector)
        i_frac = [i/denom for (i, z) in zip(i_frac_vector, self.z0) if z]
