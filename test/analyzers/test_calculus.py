import unittest
import numpy as np
import sympy as sp
from ggcas.analyzers.calculus import compute_numerical_function, compute_error, gaus_legendre_integrator, king_integrator

class TestCalculusFunctions(unittest.TestCase):

    def test_compute_numerical_function(self):
        x = sp.symbols('x')
        func = x**2
        variables = [x]
        var_data = [np.array([1, 2, 3])]
        result = compute_numerical_function(func, variables, var_data)
        expected = np.array([1, 4, 9])
        np.testing.assert_equal(result.shape, expected.shape)
        np.testing.assert_array_almost_equal(result, expected)

    def test_compute_error(self):
        x, y = sp.symbols('x y')
        func = x + y
        variables = [x, y]
        var_data = [np.array([1, 2]), np.array([3, 4])]
        var_errors = [np.array([0.1, 0.1]), np.array([0.1, 0.1])]
        result = compute_error(func, variables, var_data, var_errors)
        expected = np.array([3.1622776601683795, 4.47213595499958])
        np.testing.assert_array_almost_equal(result, expected)

    def test_gaus_legendre_integrator(self):
        def f(x):
            return x**2
        a = 0
        b = 1
        points = 20
        result = gaus_legendre_integrator(f, a, b, points)
        expected = 1/3  # Integral of x^2 from 0 to 1
        self.assertAlmostEqual(result, expected, places=5)

    def test_king_integrator(self):
        # This test assumes the Fortran90 code and required files are correctly set up.
        # It will not run correctly in this isolated environment.
        w0 = 5.0
        output = 'profile'
        result = king_integrator(w0, output)
        self.assertIsInstance(result, str)

if __name__ == '__main__':
    unittest.main()