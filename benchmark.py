import emigrate
import ionize
import cProfile

solutions = [ionize.Solution(['hepes', 'tris'],
                             [.05, .105]),
             ionize.Solution(['caproic acid', 'tris', 'fluorescein',
                              'alexa fluor 488', 'mops', 'acetic acid'],
                             [.01, .1, .01, .01, .01, .01]),
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
cProfile.run("solver.solve('examples/benchmark.hdf5', dt, tmax)")
