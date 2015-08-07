import unittest

from .Solver import Solver
from .Frame import Frame
from .FrameSeries import FrameSeries
import ionize

solutions = [ionize.Solution(['hepes', 'tris'], [.05, .105]),
             ionize.Solution(['caproic acid', 'tris',
                              'fluorescein', 'alexa fluor 488',
                              'mops', 'acetic acid'],
                             [.01, .1, .01, .01, .01, .01]),
             ionize.Solution(['hydrochloric acid', 'tris'],
                             [.04, .08]),
             ]

initialization_dict = dict(n_nodes=137,
                           lengths=[.005, .001, .02],
                           interface_length=.0005,
                           solutions=solutions,
                           current=-500.,
                           )

# #TODO:10 Add smaller unit tests.
# class TestEquilibrate(unittest.TestCase):
#     pass
#
#
# class TestFluxer(unittest.TestCase):
#     pass


class TestFrame(unittest.TestCase):
    def test_initialize(self):
        frame = Frame(initialization_dict)


class TestSolver(unittest.TestCase):
    def setUp(self):
        self.frame = Frame(initialization_dict)
        self.tmax = 20
        self.dt = 1

    def test_solve(self):
        solver = Solver(self.frame, filename='../example_1.hdf5',
                        precondition=False, flux_mode='slip')

        solver.solve(self.dt, self.tmax)

    def test_precondition(self):
        solver = Solver(self.frame, filename='../example_1.hdf5',
                        precondition=True, flux_mode='slip')

    def test_fixed_pH(self):
        self.frame.pH = 7
        solver = Solver(self.frame, filename='../example_1.hdf5',
                        precondition=False, equilibrium_mode='fixed')
        solver.solve(self.dt, self.tmax)

    def test_compact(self):
        solver = Solver(self.frame, filename='../example_1.hdf5',
                        precondition=False, flux_mode='compact')
        solver.solve(self.dt, self.tmax)



    # def test_iter(self):
    #     solver = Solver(self.frame, filename='../example_1.hdf5',
    #                     precondition=True, flux_mode='slip')
    #     tmax = 200
    #     dt = 1
    #     ode_solver = 'dopri5'
    #     for frame in solver.iterate(dt, tmax, method=ode_solver):
    #         pass

if __name__ == '__main__':
    unittest.main()
