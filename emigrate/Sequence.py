import h5py
import numpy as np
import sys
import ionize
from numbers import Number

from .__version__ import __version__
from .Frame import Frame
from .deserialize import deserialize

# Create a string data type
string_datatype = h5py.special_dtype(vlen=unicode)


class Sequence(object):

    # Public Attribute
    path = None

    # Private Attribute
    _hdf5 = None
    _compression = 'gzip'

    def __init__(self, path='default.hdf5', mode=None):
        self.path = path
        self._hdf5 = h5py.File(self.path, mode=mode)
        self.version()

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

        # If a particular dataset was serialized before storing,
        # deserialize it before using.
        for key in data.keys():
            if 'serialized' in data[key].attrs.keys():
                data[key] = [deserialize(value)
                             for value in data[key][()].tolist()]
        data.update(self._frames()[str(idx)].attrs)
        return Frame(data)

    def __iter__(self):
        def f():
            for idx in range(len(self._frames().keys())):
                yield self[idx]
        return f()

    def __setitem__(self, idx, frame):
        if not isinstance(idx, int):
            raise IndexError('Sequence index must be an integer.')

        if str(idx-1) not in self._frames().keys() and self._frames().keys():
            raise IndexError('Sequence index out of range.')

        # Create location information
        idx = str(idx)
        location = self._frames().create_group(idx)

        # Write dict items
        for key, value in frame.__dict__.items():
            if isinstance(value, np.ndarray):
                location.create_dataset(key, data=value, dtype='f4',
                                        compression=self._compression)
            elif isinstance(value, (Number, basestring)):
                location.attrs[key] = value
            else:
                try:
                    self._write_objects(key, value, location)
                except:
                    msg = "Couldn't write {}. Type: {}"
                    warnings.warn(msg.format(key, type(value)))
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

    def _write_objects(self, key, value, group):
        try:
            value = [value.serialize()]
        except:
            value = [v.serialize() for v in value]
        if sys.version_info > (3,):
            value = [v.encode('ascii') for v in value]

        data = group.create_dataset(key,
                                    data=value,
                                    shape=(len(value),),
                                    dtype=string_datatype)
        data.attrs['serialized'] = True

    def version(self):
        if 'version' not in self._hdf5.attrs.keys():
            self._hdf5.attrs['version'] = __version__
        return self._hdf5.attrs['version']

    def mode(self):
        return self._hdf5.mode
