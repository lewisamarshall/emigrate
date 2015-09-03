import unittest
import numpy as np
from scipy.special import erf
import ionize
from click.testing import CliRunner
import traceback

from .__main__ import cli
from .Solver import Solver
from .Frame import Frame
from .Sequence import Sequence
from .flux_schemes.Differentiate import Differentiate
from .deserialize import deserialize
from .equilibration_schemes.Multiroot import Multiroot

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


class TestMultiroot(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.array = np.array([[1, 1, 3, -10],
                               [1, 8, 4, -12]]
                              ).transpose()
        self.array = np.concatenate([self.array]*1000, axis=1)
        self.multiroot = Multiroot()

    def test_multiroot_noguess(self):
        self.multiroot(self.array)

    def test_multiroot_guess(self):
        guess = self.multiroot(self.array)
        self.multiroot(self.array, guess-1)


class TestDifferentiate(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.Nt = 50
        z = np.linspace(-1, 1, self.Nt)
        self.test_functions = np.array([19*erf(z*3),
                                        25*erf(z*2),
                                        (15 * (erf(z*2) +
                                         .3 * np.random.random(z.shape)) /
                                         (10*z**2+1))
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

    def test_serialize(self):
        frame = Frame(initialization_dict)
        recreated = deserialize(frame.serialize(compact=False))
        self.assertEqual(frame.serialize(compact=True),
                         recreated.serialize(compact=True))


class TestSequence(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.frame = Frame(initialization_dict)
        self.fs = Sequence(path='examples/test_sequence.hdf5', mode='w')
        self.fs.append(self.frame)

    def test_add(self):
        for time in range(5):
            self.fs.append(self.frame)

    def test_iterate(self):
        [None for frame in self.fs]


class TestSolver(unittest.TestCase):
    def setUp(self):
        self.frame = Frame(initialization_dict)
        self.tmax = 2
        self.dt = 1

    def test_slip(self):
        solver = Solver(self.frame,
                        precondition=False, flux_mode='slip')

        solver.solve('examples/test_slip.hdf5', self.dt, self.tmax)

    def test_precondition(self):
        solver = Solver(self.frame, precondition=True, flux_mode='slip')

    def test_fixed_pH(self):
        self.frame.pH = 7
        solver = Solver(self.frame,
                        precondition=False, equilibrium_mode='fixed')
        solver.solve('examples/test_fixed_pH.hdf5', self.dt, self.tmax)

    def test_compact(self):
        solver = Solver(self.frame,
                        precondition=False, flux_mode='compact')
        solver.solve('examples/test_compact.hdf5', self.dt, self.tmax)

    def test_reference_frame(self):
        solutions = [ionize.Solution(['acetic acid', 'b-alanine'],
                                     [.02, .012]),
                     ionize.Solution(['acetic acid', 'pyridine', 'aniline'],
                                     [.02, .009/2, .011/2]),
                     ionize.Solution(['acetic acid', 'sodium'],
                                     [.02, .018]),
                     ]
        system = Frame(dict(lengths=[.01, .004, .01],
                            n_nodes=100,
                            interface_length=.0005,
                            solutions=solutions,
                            current_density=500,
                            domain_mode='left',
                            ))
        solver = Solver(system,
                        precondition=True,
                        flux_mode='slip')
        solver.set_reference_frame(ionize.load_ion('acetic acid'), 'right')
        solver.solve('examples/test_reference_frame.hdf5', self.dt, self.tmax)


class TestCLI(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.runner = runner = CliRunner()
        self.runner.invoke(cli,
                                    ['construct',
                                     '-i', 'examples/constructor.json',
                                     '-o', 'examples/initial_condition.json'],
                                    obj={'sequence': None, 'frame': None})
        result = self.runner.invoke(cli,
                                    ['load', 'examples/initial_condition.json',
                                     'solve', '-t', '10.0', '-d', '1.0',
                                     '--output', 'examples/cli_test.hdf5'],
                                    obj={'sequence': None, 'frame': None})

    def setUp(self):
        self.runner = runner = CliRunner()

    def tearDown(self):
        del self.runner

    def test_construct(self):
        result = self.runner.invoke(cli,
                                    ['construct',
                                     '-i', 'examples/constructor.json',
                                     '-o', 'examples/initial_condition.json'],
                                    obj={'sequence': None, 'frame': None})
        self.assertEqual(result.exit_code,
                         0,
                         ' '.join(traceback.format_tb(result.exc_info[2]))
                         )

    def test_solve(self):
        result = self.runner.invoke(cli,
                                    ['load', 'examples/initial_condition.json',
                                     'solve', '-t', '3.0', '-d', '1.0',
                                     '--output',
                                     'examples/cli_test_solve.hdf5'],
                                    obj={'sequence': None, 'frame': None})
        self.assertEqual(result.exit_code,
                         0,
                         ' '.join(traceback.format_tb(result.exc_info[2]))
                         )

    def test_load(self):
        result = self.runner.invoke(cli,
                                    ['load', 'examples/cli_test.hdf5'],
                                    obj={'sequence': None, 'frame': None})
        self.assertEqual(result.exit_code,
                         0,
                         ' '.join(traceback.format_tb(result.exc_info[2]))
                         )

    def test_echo(self):
        result = self.runner.invoke(cli,
                                    ['echo', '-f', '5'],
                                    obj={'sequence':
                                         Sequence(path='examples/cli_test.hdf5'),
                                         'frame': None})
        self.assertEqual(result.exit_code,
                         0,
                         ' '.join(traceback.format_tb(result.exc_info[2]))
                         )

    def test_plot(self):
        result = self.runner.invoke(cli,
                                    ['plot', '-f', '1',
                                     'examples/test_plot.png'],
                                    obj={'sequence':
                                         Sequence(path='examples/cli_test.hdf5'),
                                         'frame': None})
        self.assertEqual(result.exit_code,
                         0,
                         ' '.join(traceback.format_tb(result.exc_info[2]))
                         )


if __name__ == '__main__':
    unittest.main()
