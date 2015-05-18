import h5py
import json
import numpy as np
from collections import OrderedDict
from Electrolyte import Electrolyte

# Create a string data type
string_dt = h5py.special_dtype(vlen=bytes)


class Electromigration(object):

    filename = None
    hdf5 = None
    attributes = None
    electrolytes = None
    ions = None
    compression = 'gzip'
    mode = None

    def __init__(self, ions=None, filename=None, mode='w-'):
        self.filename = filename
        self.ions = ions
        self.mode = mode

        # If there is a file, open it with mode.
        if filename is not None:
            self.hdf5 = h5py.File(filename, mode)
            self._initialize_hdf5()
        # If there is no file, store everything in memory.
        else:
            self.attributes = dict()
            self.electrolytes = OrderedDict()

        self._update_ions()

    def _initialize_hdf5(self):
        self.attributes = self.hdf5.attrs
        self.electrolytes = self.hdf5.require_group('electrolytes')
        self.hdf5.flush()

    def _update_ions(self):
        if self.ions is None:
            self._ions = self.hdf5['ions']
            self.ions = self._ions
        else:
            if self.hdf5 and 'w' in self.mode:
                self._ions = self.hdf5.create_dataset('ions',
                                                    #   data=np.array(self.ions),
                                                      (len(self.ions),),
                                                      dtype=string_dt)
                for idx, value in enumerate(self.ions):
                    self._ions[idx] = value

    def __getitem__(self, frame):
        if self.hdf5:
            frame = str(frame)
            data = dict(self.electrolytes[frame])
            data['ions'] = np.array(self.ions).tolist()
            return Electrolyte(data)
        else:
            return self.electrolytes.values()[frame]

    def __iter__(self):
        def f():
            idx = 0
            while True:
                try:
                    yield self[idx]
                    idx += 1
                except:
                    raise StopIteration
        # return (i for i in self.electrolytes.values())
        return f()

    def add_electrolyte(self, time, electrolyte):
        if self.hdf5:
            self._add_electrolyte_hdf5(time, electrolyte)
        else:
            self._add_electrolyte_dict(time, electrolyte)

    def _add_electrolyte_hdf5(self, time, electrolyte):
        # Name each electrolyte group with a string representation of a
        # consecutive integer.
        idx = str(len(self.electrolytes.keys()))
        location = self.electrolytes.create_group(idx)

        # Write to location
        location.attrs['time'] = time
        self._write_electrolyte(electrolyte, location)
        self.hdf5.flush()

    def _write_electrolyte(self, electrolyte, location):
        for key, value in electrolyte.__dict__.items():
            if key in ['concentrations', 'nodes', 'pH']:
                location.create_dataset(key, data=value,
                                        compression=self.compression,
                                        dtype='f4')
            elif key is 'ions':
                location['ions'] = self._ions
            else:
                pass
                # location.attrs[key] = value

    def _add_electrolyte_dict(self, time, electrolyte):
        if time not in self.electrolytes.keys():
            self.electrolytes[time] = electrolyte

    def _serialize_electrolytes(self):
        return [electrolyte.serialize() for electrolyte in self.electrolytes.values()]

    def serialize(self, full=False):
        serial = dict()
        serial['time'] = self.electrolytes.keys()
        serial['ions'] = [ion.name for ion in self.ions]
        serial['attributes'] = self.attributes
        serial['electrolytes'] = self._serialize_electrolytes()
        return serial

    def deserialize(self, serial):
        self.ions = serial['ions']
        self.attributes = serial['attributes']
        self._deserialize_electrolytes(serial['electrolytes'], serial['time'])

    def _deserialize_electrolytes(self, serial_electrolyte, serial_time):
        for time, electrolyte in zip(serial_time, serial_electrolyte):
            target[time] = Electrolyte().deserialize(electrolyte)

if __name__ == '__main__':
    file = '/Users/lewis/Documents/github/emigrate/test.hdf5'
    ions = [str(i) for i in range(5)]
    e = Electromigration(filename=file, mode='r')
    # print e[1].concentrations.shape
    # print e[1].nodes.shape
    # print e[2].concentrations
    # print e[3].serialize()
    for idx, electrolyte in enumerate(e):
        print electrolyte.pH
