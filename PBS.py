import emigrate
import ionize
import cProfile

solutions = [ionize.Solution(['hepes', 'tris'],
                             [.2, .505]),
             ionize.Solution(['phosphoric acid',
                              'chloride',
                              'potassium',
                              'sodium',
                              'tris',
                              ],
                             [0.012, .140, .0045, .157, .03]),
             ionize.Solution(['hydrochloric acid', 'tris'],
                             [.05, .10]),
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
dt = 1
ode_solver = 'dopri5'
profile = True
solver.solve('examples/PBS_tris.hdf5', dt, tmax)
# cProfile.run("solver.solve('examples/benchmark.hdf5', dt, tmax)")
