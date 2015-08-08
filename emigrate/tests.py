import unittest
import numpy as np
from scipy.special import erf


from .Solver import Solver
from .Frame import Frame
from .FrameSeries import FrameSeries
from .flux_schemes.Differentiate import Differentiate
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


class TestMultiroot(unittest.TestCase):
    pass


class TestDifferentiate(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.Nt = 50
        z = np.linspace(-1, 1, self.Nt)
        self.test_functions = np.array([19*erf(z*3),
                                        25*erf(z*2),
                                        (15 * (erf(z*2) +
                                         .3 * np.random.random(z.shape))
                                         / (10*z**2+1))
                                        ])

    def test_6thOrder(self):
        differ = Differentiate(self.Nt, 1, method='6th-Order')
        self.differentiate_functions(differ)

    def test_dissipative(self):
        differ = Differentiate(self.Nt, 1, method='dissipative')
        self.differentiate_functions(differ)

    def differentiate_functions(self, differ):
        differ.first_derivative(self.test_functions.T)
        differ.second_derivative(self.test_functions.T)
        differ.smooth(self.test_functions.T)


# class TestEquilibrate(unittest.TestCase):
#     pass


# class TestFluxer(unittest.TestCase):
#     pass


class TestFrame(unittest.TestCase):
    def test_initialize(self):
        frame = Frame(initialization_dict)


class TestFrameSeries(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.frame = Frame(initialization_dict)
        # TODO: make frameseries properly save ions.
        self.fs = FrameSeries(ions=[ion.name for ion in self.frame.ions],
                              filename='examples/test_frame_series.hdf5',
                              mode='w')
        self.fs.add_frame(0, self.frame)

    def test_add(self):
        for time in range(5):
            self.fs.add_frame(time, self.frame)

    def test_iterate(self):
        [None for frame in self.fs]


class TestSolver(unittest.TestCase):
    def setUp(self):
        self.frame = Frame(initialization_dict)
        self.tmax = 20
        self.dt = 1

    def test_slip(self):
        solver = Solver(self.frame, filename='examples/test_slip.hdf5',
                        precondition=False, flux_mode='slip')

        solver.solve(self.dt, self.tmax)

    def test_precondition(self):
        solver = Solver(self.frame, filename='examples/test_precondition.hdf5',
                        precondition=True, flux_mode='slip')

    def test_fixed_pH(self):
        self.frame.pH = 7
        solver = Solver(self.frame, filename='examples/test_fixed_pH.hdf5',
                        precondition=False, equilibrium_mode='fixed')
        solver.solve(self.dt, self.tmax)

    def test_compact(self):
        solver = Solver(self.frame, filename='examples/test_compact.hdf5',
                        precondition=False, flux_mode='compact')
        solver.solve(self.dt, self.tmax)

if __name__ == '__main__':
    unittest.main()
