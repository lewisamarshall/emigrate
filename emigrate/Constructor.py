import numpy as np
from scipy.special import erf
import ionize


class Constructor(object):
    domain = None
    ions = None
    concentrations = None

    def __init__(self, domain_length, nodes, interface_length,
                 solutions, V, I=None, domain_mode='centered'
                 ):
        self.domain_length = domain_length
        self.nodes = nodes
        self.interface_length = interface_length
        self.solutions = solutions
        self.V = V
        self.I = I
        self.domain_mode = domain_mode

        self.create_domain()
        self.create_ions()
        self.create_concentrations()

    def create_domain(self):
        if self.domain_mode == 'centered':
            self.domain = np.linspace(-self.domain_length/2.,
                                      self.domain_length/2.,
                                      self.nodes)

        elif self.domain_mode == 'left':
            self.domain = np.linspace(0,
                                      self.domain_length,
                                      self.nodes)

    def create_ions(self):
        self.ions = []
        for solution in self.solutions:
            self.ions.extend([ion.name for ion in solution.ions])
        self.ions = list(set(self.ions))
        self.ions = ionize.Solution(self.ions, [0]*len(self.ions)).ions

    def create_concentrations(self):
        self.concentrations = []
        for ion in self.ions:
            cs = [solution.get_concentration(ion) for solution in self.solutions]
            self.concentrations.append((erf(self.domain/self.interface_length)/2.+0.5)*
                                       (cs[1]-cs[0]) + cs[0] )
        self.concentrations = np.array(self.concentrations)

if __name__ == '__main__':
    import ionize
    from matplotlib import pyplot as plot

    tris = ionize.load_ion('tris')
    solutions = [ionize.Solution(['hydrochloric acid', tris], [.05, .1]),
                 ionize.Solution(['caproic acid', tris], [.05, .1])
                 ]
    system = Constructor(domain_length=0.1,
                         nodes=50,
                         interface_length=.05,
                         solutions=solutions,
                         V=500, I=None,
                         domain_mode='centered'
                         )
    print system.domain
    print system.ions
    print system.concentrations
    for i in range(3):
        plot.plot(system.domain, system.concentrations[i,:])
    plot.show()
