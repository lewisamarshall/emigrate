import numpy as np
from Equilibrate_Base import Equilibrate_Base
import scipy
from math import log10
# pylint: disable=W0232, E1101, W0201


class Variable_pH(Equilibrate_Base):

    """A solver that uses calculates pH at each time point."""

    def calc_equilibrium(self):
        """Calculate equilibrium."""
        self.calc_pH()
        self.calc_mobility()
        self.calc_diffusivity()
        self.calc_molar_conductivity()
        pass

    def calc_pH(self, I=0):
        """Return the pH of the object."""
        # Find the order of the polynomial. This is the maximum
        # size of the list of charge states in an ion.
        max_columns = max([max(i.z)-min(i.z)+2 for i in self.ions])
        m_ions = len(self.ions)
        nodes = self.concentrations.shape[0]

        # Set up the matrix of Ls, the multiplication
        # of acidity coefficients for each ion.
        l_matrix = np.zeros([m_ions, max_columns])

        for i in range(m_ions):
            l_matrix[i, 0:len(self.ions[i].z)+1] = self.ions[i].L(I)

        # Construct Q vector.
        Q = 1.0
        for j in range(l_matrix.shape[0]):
            Q = np.convolve(Q, l_matrix[j, :])

        # Convolve with water dissociation.
        Q = np.convolve(Q, [-self.Kw_eff(I), 0.0, 1.0])

        # Construct P matrix
        PMat = []
        for i in range(m_ions):
            z_list = self.ions[i].z0

            tmp = np.zeros([1, l_matrix.shape[1]])
            tmp[0, 0:len(z_list)] = z_list
            Mmod = l_matrix.copy()
            Mmod[i, :] = Mmod[i, :] * tmp

            Pi = 1
            for kl in range(Mmod.shape[0]):
                Pi = np.convolve(Pi, Mmod[kl, :])

            Pi = np.convolve([0.0, 1.0], Pi)  # Convolve with P2
            PMat.append(Pi)

        PMat = np.array(PMat, ndmin=2)

        # Multiply P matrix by concentrations, and sum.
        C = np.tile((self.concentrations), ((PMat.shape[1], 1))).transpose()
        P = np.sum(PMat*C, 0)

        # Construct polynomial. Change the shapes as needed.
        if len(P) < len(Q):
            P.resize(Q.shape)
        elif len(P) > len(Q):
            Q.resize(P.shape)

        poly = P + Q

        # Solve Polynomial for concentration
        # reverse order for poly function
        roo = np.roots(poly[::-1])

        roo_reduced = [r for r in roo if r.real > 0 and r.imag == 0]
        if roo_reduced:
            cH = float(roo_reduced[-1].real)
        else:
            print 'Failed to find pH.'

        # Convert to pH. Use the activity to correct the calculation.
        self.pH = -log10(cH)

    def calc_mobility(self):
        pass

    def calc_diffusivity(self):
        pass

    def calc_molar_conductivity(self):
        pass
