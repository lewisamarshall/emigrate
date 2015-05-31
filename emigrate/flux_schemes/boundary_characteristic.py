import numpy as np
from numpy.linalg import eig, inv


def boundary_characteristic(self, edge):
    if edge == 'right':
        index = -1
    elif edge == 'left':
        index = 1
    else:
        raise NotImplimentedError

    a, v, w = self._get_characteristic_matricies(index)
    # print 'w:', w

    dcdt = np.dot(a,  (self.cz[:, index] / self.xz[index]))
    dRdt = np.dot(inv(v), self.dcdt[:, index])
    if edge == 'right':
        dRdt = np.where(np.greater(0, w), dRdt, 0)
    else:
        dRdt = np.where(np.greater(w, 0), dRdt, 0)
    dcdt = np.dot(v, dRdt)
    return dcdt


def _get_characteristic_matricies(self, index):
    a = self._a_matrix(index)
    w, v = eig(a)
    return a, v, w


def _a_matrix(self, index):
    a1 = np.diag(self.mobility[:, index] * self.E[index] +
                 self.bulk_flow)
    q = self.j / self.conductivity()[index]**2 * self.concentrations[:, index] *\
        self.mobility[:, index]
    r = self.molar_conductivity[:, index]
    a2 = q[:, np.newaxis].transpose() * r
    a = a1 + a2
    return a
