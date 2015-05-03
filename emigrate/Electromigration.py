import json
from collections import OrderedDict
import numpy as np
from Electrolyte import Electrolyte


class Electromigration(object):

    properties = None
    electrolytes = None
    full_electrolytes = None
    ions = None

    def __init__(self, ions):
        self.ions = ions
        self.properties = dict()
        self.electrolytes = OrderedDict()
        self.full_electrolytes = OrderedDict()

    def add_electrolyte(self, time, electrolyte, full=False):
        if not full:
            if time not in self.electrolytes.keys():
                self.electrolytes[time] = electrolyte

        if time not in self.full_electrolytes.keys():
            self.full_electrolytes[time] = electrolyte

    def write_json(self, file, full=False):
        with open(file, 'w') as open_file:
            json.dump(self.serialize(full=False), open_file)

    def load_json(self, file):
        with open(file, 'r') as open_file:
            save_object = json.load(open_file)
        return self.deserialize(save_object)

    def _serialize_electrolytes(self, full=False):
        if full:
            target = self.full_electrolytes
        else:
            target = self.electrolytes

        return [electrolyte.serialize() for electrolyte in target.values()]

    def serialize(self, full=False):
        serial = dict()
        serial['time'] = self.electrolytes.keys()
        serial['ions'] = [ion.name for ion in self.ions]
        serial['properties'] = self.properties
        serial['electrolytes'] = self._serialize_electrolytes(full=False)
        if full:
            serial['full_electrolytes'] = self._serialize_electrolytes(full=True)
            serial['full_time'] = self.full_electrolytes.keys()
        return serial

    def deserialize(self, serial):
        self.ions = serial['ions']
        self.properties = serial['properties']
        self._deserialize_electrolytes(serial['electrolytes'], serial['time'], full=False)
        self._deserialize_electrolytes(serial['full_electrolytes'], serial['full_time'], full=True)
        return self

    def _deserialize_electrolytes(self, serial_electrolyte, serial_time, full=False):
        if full:
            target = self.full_electrolytes
        else:
            target = self.electrolytes

        for time, electrolyte in zip(serial_time, serial_electrolyte):
            target[time] = Electrolyte().deserialize(electrolyte)
