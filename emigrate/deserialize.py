import json
import numpy as np
from ionize.deserialize import object_hook as ionize_hook
import ionize
from .Frame import Frame


def deserialize(serial):
    serial = json.loads(serial, object_hook=object_hook)
    return serial


def object_hook(obj):
    if '__frame__' in obj:
        obj.pop('__frame__')
        return Frame(obj)
    elif '__ndarray__' in obj:
        return np.array(obj['data'])
    return ionize_hook(obj)
