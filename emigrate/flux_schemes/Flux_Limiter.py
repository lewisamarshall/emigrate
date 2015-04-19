"""Defines a flux limiter class."""
import numpy as np


class Flux_Limiter(object):

    """Return a flux limiting function."""

    mode = False
    limiter_dict = dict()

    def __init__(self, mode):
        """Initialize."""
        self.mode = mode
        self.limiter_dict = {"minmod": self.minmod,
                             "superbee": self.superbee,
                             "koren": self.koren
                             }

    def limit(self, solution):
        """Select a function and return a limiter."""
        if self.mode:
            r = self.get_r(solution)
            return self.limiter_dict[self.mode](r)
        else:
            return np.zeros(solution.shape)

    def get_r(self, solution):
        """Calculate r, the input parameter for limiters."""
        diff = np.diff(solution, 1)
        self.diff = diff
        r = diff[:, 1:]/diff[:, :-1]
        r = np.pad(r, ((0.,0.),(1.,1.)), 'constant', constant_values=0)
        return r

    def minmod(self, r):
        """The minmod limiter function."""
        # return np.maximum(0., np.minimum(1.,r))
        v = self.diff[:, 2:]
        w = self.diff[:, :-2]
        L = 1./2.*(v + w)*(1-np.fabs((v-w)/(np.fabs(v)-np.fabs(w))))

    def superbee(self, r):
        """The superbee limiter function."""
        return np.maximum(0, np.maximum(np.minimum(1.,2.*r), np.minimum(r,2.)))

    def koren(self,r):
        """The koren limiter function."""
        return np.maximum(0, np.minimum(np.minimum(2*r, (1+2*r)/3), 2))


if __name__ == '__main__':
    import matplotlib.pyplot as plot
    x = np.linspace(-3,3)
    y = Flux_Limiter('superbee').koren(x)
    plot.plot(x,y)
    plot.show()
