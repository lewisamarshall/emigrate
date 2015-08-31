import ionize
import simplejson as json

solutions = [ionize.Solution(['hepes', 'tris'],
                             [.05, .105]),
             ionize.Solution(['caproic acid', 'tris', 'fluorescein',
                              'alexa fluor 488', 'mops', 'acetic acid'],
                             [.01, .1, .01, .01, .01, .01]),
             ionize.Solution(['hydrochloric acid', 'tris'],
                             [.04, .08]),
             ]

constructor = {'n_nodes': 137,
               'lengths': [.005, .001, .02],
               'interface_length': .0005,
               'solutions': solutions,
               'current': -500.,
               }


def _encode(obj):
    if isinstance(obj, ionize.Ion):
        ion = obj.serialize(nested=True)
        return ion
    elif isinstance(obj, ionize.Solution):
        solution = obj.serialize(nested=True)
        return solution
    elif isinstance(obj, np.ndarray):
        return {'__ndarray__': True, 'data': obj.tolist()}
    elif isinstance(obj, h5py._hl.dataset.Dataset):
        return {'__ndarray__': True, 'data': obj[()].tolist()}
    return json.JSONEncoder().default(obj)

with open('constructor.json', 'w') as fileobj:
    json.dump(constructor, fileobj, default=_encode, sort_keys=True,
              indent=4, separators=(',', ': '))
