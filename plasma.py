import emigrate
import ionize
import cProfile

solutions = [ionize.Solution(['hepes', 'tris'],
                             [.2, .505]),
             ionize.Solution(['carbonic acid',
                              'chloride',
                              'potassium',
                              'sodium',
                              'tris',
                              'alexa fluor 488',
                              ],
                             [0.032, .106, .006, .132, .03, .0001]),
             ionize.Solution(['hydrochloric acid', 'tris'],
                             [.05, .10]),
             ]

system = emigrate.Frame(dict(
                         n_nodes=137,
                         lengths=[.005, .02, .04],
                         interface_length=.0005,
                         solutions=solutions,
                         current=-500.,
                         ))

solver = emigrate.Solver(system,
                         precondition=True,
                         flux_mode='slip')
tmax = 400
dt = 1
ode_solver = 'dopri5'
profile = True
solver.solve('examples/plasma.hdf5', dt, tmax)
# cProfile.run("solver.solve('examples/benchmark.hdf5', dt, tmax)")
