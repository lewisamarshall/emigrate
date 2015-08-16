import json
import numpy as np
# from ionize.deserialize import object_hook as ionize_hook
import ionize
from .Frame import Frame
from .FrameSeries import FrameSeries

def deserialize(serial):
    serial = json.loads(serial, object_hook=object_hook)
    return serial

def object_hook(obj):
    if '__frame__' in obj:
        obj.pop('__frame__')
        return Frame(obj)
    elif '__ndarray__' in obj:
        return np.array(obj['data'])
    elif '__ion__' in obj:
        obj.pop('__ion__')
        obj.pop('type')
        return ionize.Ion(**obj)
    elif '__solution__' in obj:
        obj.pop('__solution__')
        obj.pop('type')
        return ionize.Solution(**obj)
    elif 'type' in obj:
        t = obj.pop('type')
        if 'ion' in t:
            return ionize.Ion(**obj)
    return obj
