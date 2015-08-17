import h5py
import json
import numpy as np
import sys
from collections import OrderedDict
from Frame import Frame

# Create a string data type
string_dt = h5py.special_dtype(vlen=unicode)


class FrameSeries(object):

    filename = None
    hdf5 = None
    attributes = None
    frames = None
    ions = None
    compression = 'gzip'
    mode = None

    def __init__(self, ions=None, filename='default.hdf5', mode='w-'):
        self.filename = filename
        self.ions = ions
        self.mode = mode

        # If there is a file, open it with mode.
        self.hdf5 = h5py.File(filename, mode)
        self._initialize_hdf5()
        self._update_ions()

    def _initialize_hdf5(self):
        self.attributes = self.hdf5.attrs
        self.frames = self.hdf5.require_group('frames')
        self.hdf5.flush()

    def _update_ions(self):
        if self.ions is None:
            self._ions = self.hdf5['ions']
            self.ions = self._ions
        else:
            if self.hdf5 and 'w' in self.mode:
                self._ions = self.hdf5.create_dataset('ions',
                                                      (len(self.ions),),
                                                      dtype=string_dt)
                for idx, value in enumerate(self.ions):
                    if sys.version_info < (3,):
                        self._ions[idx] = value
                    else:
                        self._ions[idx] = value.encode('ascii')

    def __getitem__(self, frame):
        frame = str(frame)
        data = dict(self.frames[frame])
        data['ions'] = np.array(self.ions).tolist()
        return Frame(data)

    def __iter__(self):
        def f():
            idx = 0
            while True:
                try:
                    yield self[idx]
                    idx += 1
                except:
                    raise StopIteration
        return f()

    def add_frame(self, time, frame):
        # Name each electrolyte group with a string representation of a
        # consecutive integer.
        idx = str(len(self.frames.keys()))
        location = self.frames.create_group(idx)

        # Write to location
        location.attrs['time'] = time
        self._write_frame(frame, location)
        self.hdf5.flush()

    def _write_frame(self, frame, location):
        for key, value in frame.__dict__.items():
            if key in ['concentrations', 'nodes', 'pH', 'field']:
                try:
                    location.create_dataset(key, data=value,
                                            compression=self.compression,
                                            dtype='f4')
                except TypeError:
                    pass
            elif key is 'ions':
                location['ions'] = self._ions
            else:
                pass


if __name__ == '__main__':
    file = '/Users/lewis/Documents/github/emigrate/test.hdf5'
    ions = [str(i) for i in range(5)]
    e = FrameSeries(filename=file, mode='r')
    # print e[1].concentrations.shape
    # print e[1].nodes.shape
    # print e[2].concentrations
    # print e[3].serialize()
    for idx, electrolyte in enumerate(e):
        print electrolyte.pH
