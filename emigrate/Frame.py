import numpy as np
from scipy.special import erf
import ionize
import warnings


class Frame(object):

    """Represent an electrophoresis system."""

    # Required Properties
    nodes = None
    ions = []
    concentrations = None

    # Optional Properties
    pH = None
    ionic_strength = None

    # Default Properties
    voltage = 0
    current_density = 0
    current = 0.
    area = 1.
    pressure = 0
    bulk_flow = 0

    def __init__(self, constructor):
        """Initialize a Frame object."""

        # Try loading from file first.
        if isinstance(constructor, basestring):
            self.open_file(constructor)

        # Next look for a construction dictionary
        elif 'solutions' in constructor.keys():
            self.construct(constructor)

        # Finally try to deserialize the object directly
        else:
            self.deserialize(constructor)

        self._resolve_current()

    def deserialize(self, constructor):
        """Use a dictionary to update electrolyte properties."""
        # Initialize optional properties
        self.ions = constructor['ions']
        for key, value in constructor.items():
            if key in self.__class__.__dict__.keys():
                if key in ['nodes', 'pH']:
                    value = np.array(value)
                elif key == 'concentrations':
                    if isinstance(value, dict):
                        value = np.array([value[i] for i in self.ions])
                    else:
                        value = np.array(value)
                # print key, value
                self.__dict__[key] = value
            else:
                warnings.warn('Unused initialization key: {}'.format(key))

    def serialize(self):
        serial = self.__dict__.copy()
        try:
            serial['ions'] = [ion.name for ion in self.ions]
        except:
            serial['ions'] = self.ions
        serial['nodes'] = self.nodes.tolist()
        serial['concentrations'] = {i: c for i, c in
                                    zip(serial['ions'],
                                        self.concentrations.tolist()
                                        )
                                    }
        serial['concentrations']['x'] = serial['nodes']
        try:
            serial['pH'] = serial['pH'].tolist()
        except:
            pass
        serial['properties'] = {'pH': serial['pH']}
        serial['properties']['x'] = serial['nodes']

        return serial

    def construct(self, constructor_input):
        """Construct electrophoretic system based on a set of solutions."""

        # Set up properties
        constructor = {'solutions': None,
                       'lengths': None,
                       'n_nodes': 100,
                       'interface_length': 1e-4,
                       'voltage': 0,
                       'current_density': 0,
                       'current': 0,
                       'domain_mode': 'left',
                       'bulk_flow': 0.,
                       'area': None,
                       }

        constructor.update(constructor_input)

        self._create_ions(constructor)
        self._create_domain(constructor)
        self._create_concentrations(constructor)
        for feature in ['current', 'current_density', 'voltage']:
            self.__dict__[feature] = constructor[feature]

    def _create_domain(self, constructor):
        """Initially place grid points in the domain."""
        constructor['domain_length'] = sum(constructor['lengths'])
        if constructor['domain_mode'] == 'centered':
            left, right = (-constructor['domain_length']/2.,
                           constructor['domain_length']/2.)

        elif constructor['domain_mode'] == 'left':
            left, right = (0, constructor['domain_length'])

        elif constructor['domain_mode'] == 'right':
            left, right = (-constructor['domain_length'], 0)

        elif constructor['domain_mode'] == 'detector':
            left, right = (-constructor.detector_location,
                           (constructor.domain_length -
                            constructor.detector_location))
        else:
            raise NotImplementedError

        self.nodes = np.linspace(left, right, constructor['n_nodes'])

    def _create_ions(self, constructor):
        self.ions = []
        for solution in constructor['solutions']:
            self.ions.extend([ion.name for ion in solution.ions])
        self.ions = list(set(self.ions))

        # Replace strings with ion objects.
        self.ions = ionize.Solution(self.ions, [0.1]*len(self.ions)).ions

    def _create_concentrations(self, constructor):
        self.concentrations = []

        for ion in self.ions:
            ion_concentration = np.zeros(self.nodes.shape)
            cs = [solution.get_concentration(ion)
                  for solution in constructor['solutions']]

            for idx in range(len(cs)):
                if idx == 0:
                    ion_concentration += cs[0]
                    left_side, right_side = (self.nodes[0],
                                             self.nodes[0] +
                                             constructor['lengths'][idx])
                else:
                    left_side, right_side = (right_side,
                                             right_side +
                                             constructor['lengths'][idx])
                    ion_concentration += ((cs[idx]-cs[idx-1]) *
                                          (erf((self.nodes-left_side)
                                           / constructor['interface_length'])
                                           / 2. + .5))

            self.concentrations.append(ion_concentration)
        self.concentrations = np.array(self.concentrations)

    def _resolve_current(self):
        if self.current:
            cd = self.current/self.area
            if self.current_density and not self.current_density == cd:
                warnings.warn('Current and current density do not match. '
                              + 'Correcting current density.')
            self.current_density = cd
        elif self.current_density:
            self.current = self.area * self.current_density

    def open_file(self, filename):
        with open(filename, 'r') as file:
            raise NotImplementedError


if __name__ == '__main__':
    my_solutions = [ionize.Solution(['hydrochloric acid', 'tris'], [.05, .1]),
                    ionize.Solution(['caproic acid', 'tris'], [.05, .1])
                    ]
    system = Electrolyte({'lengths': [0.01]*2,
                         'solutions': my_solutions})
    print system.concentrations
    print Electrolyte(system.serialize()).serialize()['concentrations'].keys()
    print Electrolyte(system.serialize()).concentrations
    # # print system.domain
    # print system.ions
    # print system.concentrations
    # for i in range(3):
    #     plot.plot(system.domain, system.concentrations[i, :])
    # plot.show()
