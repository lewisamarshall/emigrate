import numpy as np
import h5py
from scipy.special import erf
import ionize
import warnings
try:
    import simplejson as json
except:
    import json
import ionize


class Frame(object):

    """Represent an electrophoresis system."""

    # Core Properties
    nodes = None
    ions = None
    concentrations = None
    time = None

    # Solution Properties
    pH = None
    cH = None
    ionic_strength = None

    # Default Properties
    voltage = 0
    field = None
    current_density = 0
    current = 0.
    area = 1.
    pressure = 0
    bulk_flow = 0

    def __init__(self, constructor):
        """Initialize a Frame object."""

        # Next look for a construction dictionary
        if 'solutions' in constructor.keys():
            self.construct(constructor)

        else:
            self.__dict__ = constructor

        self._resolve_current()

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
                                          (erf((self.nodes-left_side) /
                                           constructor['interface_length']) /
                                           2. + .5))

            self.concentrations.append(ion_concentration)
        self.concentrations = np.array(self.concentrations)

    def _resolve_current(self):
        if self.current:
            cd = self.current/self.area
            if self.current_density and not self.current_density == cd:
                warnings.warn('Current and current density do not match. '
                              'Correcting current density.')
            self.current_density = cd
        elif self.current_density:
            self.current = self.area * self.current_density

    # I/O
    def serialize(self, compact=True):
        serial = {'__frame__': True}
        serial.update(self.__dict__)

        if compact:
            sort_keys, indent, separators = True, None, (',', ':')
        else:
            sort_keys, indent, separators = True, 4, (', ', ': ')

        return json.dumps(serial, default=self._encode, sort_keys=sort_keys,
                          indent=indent, separators=separators)

    def _encode(self, obj):
        if isinstance(obj, ionize.Ion):
            ion = obj.serialize(nested=True)
            ion.update({'__ion__': True})
            return ion
        elif isinstance(obj, np.ndarray):
            return {'__ndarray__': True, 'data': obj.tolist()}
        elif isinstance(obj, h5py._hl.dataset.Dataset):
            return {'__ndarray__': True, 'data': obj[()].tolist()}
        # Let the base class default method raise the TypeError
        return json.JSONEncoder().default(obj)


if __name__ == '__main__':
    my_solutions = [ionize.Solution(['hydrochloric acid', 'tris'], [.05, .1]),
                    ionize.Solution(['caproic acid', 'tris'], [.05, .1])
                    ]
    system = Frame({'lengths': [0.01]*2,
                    'solutions': my_solutions})
    print system.concentrations
    print system
    print system.serialize(compact=False)
