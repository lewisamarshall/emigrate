import h5py
import numpy as np
import sys
import contextlib
from .Frame import Frame

# Create a string data type
string_datatype = h5py.special_dtype(vlen=unicode)


class FrameSeries(object):

    filename = None
    hdf5 = None
    ions = None

    _compression = 'gzip'
    _attributes = None

    def __init__(self, ions=None, filename='default.hdf5'):
        self.filename = filename
        self.ions = ions

        # If there is a file, open it with mode.
        self.hdf5 = h5py.File(filename)
        self.mode('r')
        self._update_ions()

    def _frames(self):
        return self.hdf5.require_group('frames')

    def _update_ions(self):
        if self.ions is None:
            self._ions = self.hdf5['ions']
            self.ions = self._ions
        else:
            with self.mode('w'):
                self._ions = self.hdf5.create_dataset('ions',
                                                      (len(self.ions),),
                                                      dtype=string_datatype)
                for idx, value in enumerate(self.ions):
                    if sys.version_info < (3,):
                        self._ions[idx] = value
                    else:
                        self._ions[idx] = value.encode('ascii')

    def __getitem__(self, idx):
        idx = str(idx)
        data = dict(self._frames()[idx])
        data['ions'] = np.array(self.ions).tolist()
        return Frame(data)

    def __setitem__(self, idx, frame):
        idx = str(idx)
        location = self._frames().create_group(idx)
        # location.attrs['time'] = time
        self._write_frame(frame, location)
        self.hdf5.flush()


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

    def append(self, time, frame):
        idx = len(self._frames().keys())
        self[idx]=frame

    def _write_frame(self, frame, location):
        for key, value in frame.__dict__.items():
            if key in ['concentrations', 'nodes', 'pH', 'field']:
                try:
                    location.create_dataset(key, data=value,
                                            compression=self._compression,
                                            dtype='f4')
                except TypeError:
                    pass
            elif key is 'ions':
                location['ions'] = self.hdf5['ions']
            else:
                pass

    def mode(self, mode=None):
        last_mode = self.hdf5.mode
        if mode is None:
            return last_mode
        elif last_mode != mode:
                self.hdf5.close()
                self.hdf5 = h5py.File(self.filename, mode)
        return self._mode_context(last_mode)

    @contextlib.contextmanager
    def _mode_context(self, last_mode):
        yield
        self.mode(last_mode)


if __name__ == '__main__':
    filename = '/Users/lewis/Documents/github/emigrate/test.hdf5'
    ions = [str(i) for i in range(5)]
    e = FrameSeries(filename=filename, mode='r')
    # print e[1].concentrations.shape
    # print e[1].nodes.shape
    # print e[2].concentrations
    # print e[3].serialize()
    for idx, electrolyte in enumerate(e):
        print electrolyte.pH
