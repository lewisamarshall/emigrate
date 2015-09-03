import h5py
import sys
import ionize

from .Frame import Frame
from .deserialize import deserialize

# Create a string data type
string_datatype = h5py.special_dtype(vlen=unicode)


class Sequence(object):

    # Public Attribute
    path = None

    # Private Attribute
    _ions = None
    _hdf5 = None
    _compression = 'gzip'

    def __init__(self, path='default.hdf5', mode=None):
        self.path = path
        self._hdf5 = h5py.File(self.path, mode=mode)
        if 'ions' in self._hdf5.keys():
            self._ions = self._hdf5['ions']

    def __enter__(self):
        return self

    def __exit__(self, eType, eValue, eTrace):
        self.close()

    def __getitem__(self, idx):
        assert isinstance(idx, int), "Index must be an integer."
        idx = str(idx)
        assert idx in self._frames().keys(), ("Frame {} doesn't exist.").format(idx)
        data = dict(self._frames()[idx])
        data['ions'] = [deserialize(ion) for ion in self._ions[()].tolist()]
        return Frame(data)

    def __iter__(self):
        def f():
            idx = 0
            for idx in range(len(self._frames().keys())):
                yield self[idx]
            raise StopIteration

        return f()

    def __setitem__(self, idx, frame):
        assert isinstance(idx, int), "Index must be an integer."
        if len(self._frames().keys()) == 0:
            self._set_ions(frame.ions)

        # Create location information
        idx = str(idx)
        location = self._frames().create_group(idx)

        # Write dict items
        for key, value in frame.__dict__.items():
            if key in ['concentrations', 'nodes', 'pH', 'field']:
                try:
                    location.create_dataset(key, data=value,
                                            compression=self._compression,
                                            dtype='f4')
                except TypeError:
                    pass
            elif key is 'ions':
                location['ions'] = self._ions
            else:
                pass

        self._flush()

    def __len__(self):
        return len(self._frames().keys())

    def append(self, frame):
        idx = len(self)
        self[idx] = frame

    def close(self):
        self._hdf5.close()

    def _flush(self):
        self._hdf5.flush()

    def _frames(self):
        return self._hdf5.require_group('frames')

    def _set_ions(self, ions):
        self._ions = self._hdf5.create_dataset('ions',
                                               shape=(len(ions),),
                                               dtype=string_datatype)
        for idx, ion in enumerate(ions):
            if sys.version_info < (3,):
                serial = ion.serialize()
            else:
                serial = ion.serialize().encode('ascii')

            self._ions[idx] = serial

        self._flush()


if __name__ == '__main__':
    path = '/Users/lewis/Documents/github/emigrate/test.hdf5'
    ions = [str(i) for i in range(5)]
    e = Sequence(path=path, mode='r')
    # print e[1].concentrations.shape
    # print e[1].nodes.shape
    # print e[2].concentrations
    # print e[3].serialize()
    for idx, electrolyte in enumerate(e):
        print electrolyte.pH
