import emigrate
import ionize
import cProfile

solutions = [ionize.Solution(['hepes', 'tris'],
                             [.1, .2]),
             ionize.Solution(['citric acid',
                              'chloride',
                              'sodium',
                              'tris',
                              'alexa fluor 488',
                              'carbonic acid',
                              ],
                             [0.0005, .0285, .0315, .015, .0001, .01]),
             ionize.Solution(['hydrochloric acid', 'tris'],
                             [.085, .17]),
             ]

system = emigrate.Frame(dict(
                         n_nodes=137,
                         lengths=[.005, .02, .04],
                         interface_length=.0005,
                         solutions=solutions,
                         current=-400.,
                         ))

solver = emigrate.Solver(system,
                         precondition=False,
                         flux_mode='slip')
tmax = 600
dt = .2
ode_solver = 'dopri5'
profile = True
solver.solve('examples/TD_CO2.hdf5', dt, tmax)
# cProfile.run("solver.solve('examples/benchmark.hdf5', dt, tmax)")
