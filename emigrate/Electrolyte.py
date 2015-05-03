"""Contains the emigrate Constructor class."""
import numpy as np
from scipy.special import erf
import ionize


class Electrolyte(object):

    """Represent an electrophoresis system.

    Args:
        nodes=None, ions=None, concentrations=np.array(),
        pH=None, ionic_strength=None, voltage=None, current_density=None

    Attributes:
        nodes
        concentrations
        pH
        ionic_strength
        voltage
        current_density

    Raises:
        None

    Example:
        To to initialize a system, call as:

        >>> system = emigrate.Electrolyte(nodes=None,
                                          ions=None,
                                          concentrations=np.array(),
                                          pH=None,
                                          ionic_strength=None,
                                          voltage=None,
                                          current_density=None)

        To initialize a system with smooth boundaries between different
        concentration regions, use the Electrolyte.construct method.

        >>> system = emigrate.Electrolyte.construct(domain_length=0.01,
                         nodes=1005,
                         interface_length=.0005,
                         solutions=solutions,
                         V= -25., I=None,
                         domain_mode='centered'
                         )
    """

    nodes = None
    ions = []
    concentrations = None
    pH = None
    ionic_strength = None
    voltage = 0
    current_density = 0
    # current = None
    # area = None
    # pressure = 0
    # bulk_flow = 0

    def __init__(self, nodes=None, ions=None, concentrations=None,
                 pH=None, ionic_strength=None,
                 voltage=None, current_density=None, bulk_flow=0.):
        """Initialize a Electrolyte object."""
        self.nodes = nodes
        self.ions = ions
        self.concentrations = concentrations
        self.pH = pH
        self.ionic_strength = ionic_strength
        self.voltage = voltage
        self.u = bulk_flow
        self.current_density = current_density

    def construct(self, solutions, lengths, n_nodes=100,
                  interface_length=1e-4, voltage=0, current_density=0,
                  domain_mode='centered', bulk_flow = 0.):
        """Construct electrophoretic system based on a set of solutions."""
        self.voltage = voltage
        self.current_density = current_density
        self.create_ions(solutions)
        self.u = bulk_flow

        domain_length = sum(lengths)
        self.create_domain(n_nodes, domain_length, domain_mode)
        self.create_concentrations(solutions, lengths, interface_length)
        return self

    def create_domain(self, n_nodes=100, domain_length=1e-2,
                      domain_mode='centered', detector_location=None):
        """Initially place grid points in the domain."""
        if domain_mode == 'centered':
            self.nodes = np.linspace(-domain_length/2.,
                                     domain_length/2.,
                                     n_nodes)

        elif domain_mode == 'left':
            self.nodes = np.linspace(0,
                                     domain_length,
                                     n_nodes)
        elif domain_mode == 'right':
            self.nodes = np.linspace(-domain_length, 0, n_nodes)

        elif domain_mode == 'detector':
            self.nodes = np.linspace(-detector_location,
                                     domain_length-detector_location,
                                     n_nodes)
        else:
            raise NotImplementedError

    def create_ions(self, solutions):
        self.ions = []
        for solution in solutions:
            self.ions.extend([ion.name for ion in solution.ions])
        self.ions = list(set(self.ions))

        # Replace strings with ion objects.
        self.ions = ionize.Solution(self.ions, [0]*len(self.ions)).ions

    def create_concentrations(self, solutions, lengths, interface_length):
        self.concentrations = []

        for ion in self.ions:
            ion_concentration = np.zeros(self.nodes.shape)
            cs = [solution.get_concentration(ion) for solution in solutions]

            for idx in range(len(cs)):
                if idx == 0:
                    ion_concentration += cs[0]
                    left_side, right_side = (self.nodes[0],
                                             self.nodes[0] + lengths[idx])
                else:
                    left_side, right_side = (right_side,
                                             right_side + lengths[idx])
                    ion_concentration += ((cs[idx]-cs[idx-1]) *
                                          (erf((self.nodes-left_side)
                                           / interface_length) / 2. + .5))

            self.concentrations.append(ion_concentration)
        self.concentrations = np.array(self.concentrations)

    def serialize(self):
        serial = dict()
        serial['ions'] = [ion.name for ion in self.ions]
        serial['nodes'] = self.nodes.tolist()
        serial['concentrations'] = self.concentrations.tolist()
        serial['pH'] = self.pH.tolist()
        return serial

    def deserialize(self, serial):
        self.ions = serial['ions']
        self.nodes = np.array(serial['nodes'])
        self.concentrations = np.array(serial['concentrations'])
        return self

if __name__ == '__main__':
    from matplotlib import pyplot as plot

    tris = ionize.load_ion('tris')
    my_solutions = [ionize.Solution(['hydrochloric acid', tris], [.05, .1]),
                    ionize.Solution(['caproic acid', tris], [.05, .1])
                    ]
    system = Constructor(domain_length=0.1,
                         nodes=50,
                         interface_length=.05,
                         solutions=my_solutions,
                         V=500, I=None,
                         domain_mode='centered'
                         )
    print system.domain
    print system.ions
    print system.concentrations
    for i in range(3):
        plot.plot(system.domain, system.concentrations[i, :])
    plot.show()
