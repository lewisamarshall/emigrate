import json
from collections import OrderedDict
import numpy as np
from Electrolyte import Electrolyte


class Electromigration(obj):

    properties = None
    electrolytes = None
    full_electrolytes = None
    ions = None

    def __init__(self, ions):
        ions = ions
        properties = dict()
        electrolytes = OrderedDict()
        full_electrolytes = OrderedDict()

    def add_electrolyte(self, time, electrolyte, full=False):
        if not full:
            if time not in self.electrolytes.keys():
                self.electrolytes[time] = electrolyte

        if time not in self.full_electrolytes.keys():
            self.full_electrolytes[time] = electrolyte

    def write_json(self, file):
        save_object = dict()
        save_object['ions'] = [ion.name for ion in self.ions]
        save_object['properties'] = self.properties
        save_object['electrolytes'] = \
            self._serialize_electrolyte_dict(full=False)
        save_object['full_electrolytes'] = \
            self._serialize_electrolyte_dict(full=True)
        with open(file, 'w') as open_file:
            json.dump(open_file, save_object)

    def load_json(self, file):
        with open(file, 'r') as open_file:
            save_object = json.load(open_file)
        self.peroperties = save_object.properties
        self.ions = save_object.ions
        for key, nodes, concentrations in save_object['electrolytes']:
            self.electrolytes[key] = \
                Electrolyte(nodes=nodes,
                            ions=self.ions,
                            concentrations=np.array(concentrations)
                            )

    def _serialize_electrolyte_dict(self, full=False):
        if full:
            target = self.full_electrolytes
        else:
            target = self.electrolytes

        return [(time,
                 electrolyte.nodes.tolist(),
                 electrolyte.concentrations.tolist()
                 )
                for time, electrolyte in
                self.electrolytes.items()]
