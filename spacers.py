import emigrate
import ionize
import cProfile

solutions = [ionize.Solution(['hepes', 'tris',],
                             [.2, .505]),
             ionize.Solution(['tris',
                              'mops',
                              'caproic acid',
                              'ascorbic acid',
                              'carbonic acid',
                              ],
                             [.08, .01, .01, .01, .01]),
             ionize.Solution(['hydrochloric acid', 'tris'],
                             [.04, .08]),
             ]

system = emigrate.Frame(dict(
                         n_nodes=137,
                         lengths=[.005, .001, .02],
                         interface_length=.0005,
                         solutions=solutions,
                         current=-500.,
                         ))

solver = emigrate.Solver(system,
                         precondition=True,
                         flux_mode='slip')
tmax = 200
dt = 1
ode_solver = 'dopri5'
profile = True
solver.solve('examples/spacers.hdf5', dt, tmax)
