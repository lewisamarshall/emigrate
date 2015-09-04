import h5py
import sys
import ionize

from .__version__ import __version__
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
        self.version()
        if 'ions' in self._hdf5.keys():
            self._ions = self._hdf5['ions']

    def __enter__(self):
        return self

    def __exit__(self, eType, eValue, eTrace):
        self.close()

    def __getitem__(self, idx):
        if not isinstance(idx, int):
            raise IndexError('Sequence index must be an integer.')
        if not str(idx) in self._frames().keys():
            raise IndexError('Sequence index out of range.')

        data = dict(self._frames()[str(idx)])
        data['ions'] = [deserialize(ion)
                        for ion in self._ions[()].tolist()]
        return Frame(data)

    def __iter__(self):
        def f():
            idx = 0
            for idx in range(len(self._frames().keys())):
                yield self[idx]
            raise StopIteration

        return f()

    def __setitem__(self, idx, frame):
        if not isinstance(idx, int):
            raise IndexError('Sequence index must be an integer.')

        if str(idx-1) not in self._frames().keys() and self._frames().keys():
            raise IndexError('Sequence index out of range.')

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

    def __repr__(self):
        return "Sequence(path='{}', mode='{}')".format(self.path, self.mode())

    def __str__(self):
        lines = []
        lines.append('Sequence')
        lines.append('-----------------')
        lines.append("Path:    '{}'".format(self.path))
        lines.append("Mode:    '{}'".format(self.mode()))
        lines.append('Length:  {}'.format(len(self)))
        lines.append('Version: {}'.format(self.version()))
        return '\n'.join(lines)

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

    def version(self):
        if 'version' not in self._hdf5.attrs.keys():
            self._hdf5.attrs['version'] = __version__
        return self._hdf5.attrs['version']

    def mode(self):
        return self._hdf5.mode


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
