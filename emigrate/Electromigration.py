import json
from collections import OrderedDict
import numpy as np
from Electrolyte import Electrolyte
import h5py


class Electromigration(object):

    properties = None
    electrolytes = None
    full_electrolytes = None
    ions = None
    hdf5 = None

    def __init__(self, ions, hdf5=False):
        self.ions = ions
        if hdf5:
            self.create_hdf5_structure(hdf5)
        else:
            self.properties = dict()
            self.electrolytes = OrderedDict()
            self.full_electrolytes = OrderedDict()

    def add_electrolyte(self, time, electrolyte, full=False):
        if self.hdf5:
            self._add_electrolyte_hdf5(time, electrolyte, full)
        else:
            self._add_electrolyte_json(time, electrolyte, full)

    # HDF5 functions
    #####################################################################
    def create_hdf5_structure(self, filename):
        self.hdf5 = h5py.File(filename, 'w')
        self.electrolytes = self.hdf5.create_group('electrolytes')
        self.full_electrolytes = self.hdf5.create_group('full_electrolytes')
        self.hdf5.flush()

    def _add_electrolyte_hdf5(self, time, electrolyte, full=False):
        if full:
            target = self.full_electrolytes
            return None
        else:
            target = self.electrolytes

        # Name each electrolyte group with a string representation of a
        # consecutive integer.
        i = len(target.keys())+1
        electrolyte_location = target.create_group(str(i))

        # Set attributes
        electrolyte_location.attrs['time'] = time
        electrolyte_location.attrs['current density'] = \
            electrolyte.current_density
        electrolyte_location.attrs['voltage'] = electrolyte.voltage
        electrolyte.write_to_hdf5(electrolyte_location)

        self.hdf5.flush()

    def __exit__(self):
        # Close the hdf5 item on exit if it exists.
        if self.hdf5:
            self.hdf5.close()

    # JSON functions:
    #####################################################################
    def _add_electrolyte_json(self, time, electrolyte, full=False):
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
