import unittest
import numpy as np
import sympy as sp
from astropy import units as u
from ggcas.functions import *

class TestFunctions(unittest.TestCase):

    def test_angular_separation(self):
        result = angular_separation(ra0=10*u.deg, dec0=20*u.deg)
        self.assertIn('f', result)
        self.assertIn('vars', result)
        self.assertEqual(len(result['vars']), 2)

    def test_los_distance(self):
        result = los_distance()
        self.assertIn('f', result)
        self.assertIn('vars', result)
        self.assertEqual(result['vars'], sp.symbols('omega'))

    def test_radial_distance_2d(self):
        result = radial_distance_2d(analytical_w=True, ra0=10, dec0=20)
        self.assertIn('f', result)
        self.assertIn('vars', result)
        self.assertGreater(len(result['vars']), 1)

    def test_radial_distance_3d(self):
        result = radial_distance_3d(gc_distance=10*u.pc) #, analytical_r2d=True, analytical_w=True)
        self.assertIn('f', result)
        self.assertIn('vars', result)
        self.assertIn('d', result)
        self.assertGreater(len(result['vars']), 1)

    def test_total_velocity(self):
        result = total_velocity()
        self.assertIn('f', result)
        self.assertIn('vars', result)
        self.assertEqual(len(result['vars']), 2)

    def test_effective_potential(self):
        result = effective_potential(shell=True)
        self.assertIn('f', result)
        self.assertIn('vars', result)
        self.assertIn('x', result['vars'])
        self.assertIn('const', result['vars'])

if __name__ == '__main__':
    unittest.main()