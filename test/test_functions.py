import unittest
import numpy as np
from astropy import units as u
from ggcas.functions import (
    EffectivePotential,
    AngularSeparation, 
    RadialDistance3D, 
    RadialDistance2D, 
    TotalVelocity, 
    LosDistance
)

class TestAngularSeparation(unittest.TestCase):
    def setUp(self):
        self.ra0 = 10.0 * u.deg
        self.dec0 = 20.0 * u.deg
        self.angular_separation = AngularSeparation(self.ra0, self.dec0)

    def test_compute(self):
        data = [np.array([15.0]), np.array([25.0])]
        result = self.angular_separation.compute(data)
        self.assertIsInstance(result, AngularSeparation)
        np.testing.assert_array_almost_equal(result.values, np.array([6.80561008])) 

    def test_compute_error(self):
        data = [np.array([15.0]), np.array([25.0])]
        errors = [np.array([0.1]), np.array([0.1])]
        result = self.angular_separation.compute_error(data, errors)
        self.assertIsInstance(result, AngularSeparation)
        np.testing.assert_array_almost_equal(result.errors, np.array([0.09563928]))

class TestLosDistance(unittest.TestCase):
    def setUp(self):
        self.los_distance = LosDistance()

    def test_compute(self):
        data = [np.array([0.1])]
        result = self.los_distance.compute(data)
        self.assertIsInstance(result, LosDistance)
        np.testing.assert_array_almost_equal(result.values, np.array([10.]))

    def test_compute_error(self):
        data = [np.array([0.1])]
        errors = [np.array([0.01])]
        result = self.los_distance.compute_error(data, errors)
        self.assertIsInstance(result, LosDistance)
        np.testing.assert_array_almost_equal(result.errors, np.array([1.]))

class TestRadialDistance2D(unittest.TestCase):
    def setUp(self):
        self.gc_distance = 1000 * u.pc
        self.radial_distance_2d = RadialDistance2D(self.gc_distance)

    def test_compute(self):
        data = [np.array([0.1])]
        result = self.radial_distance_2d.compute(data)
        self.assertIsInstance(result, RadialDistance2D)
        np.testing.assert_array_almost_equal(result.values, np.array([100.33467209]))

    def test_compute_error(self):
        data = [np.array([0.1])]
        errors = [np.array([0.01])]
        result = self.radial_distance_2d.compute_error(data, errors)
        self.assertIsInstance(result, RadialDistance2D)
        np.testing.assert_array_almost_equal(result.errors, np.array([10.10067046]))

class TestRadialDistance3D(unittest.TestCase):
    def setUp(self):
        self.gc_distance = 1000 * u.pc
        self.radial_distance_3d = RadialDistance3D(gc_distance=self.gc_distance)

    def test_compute(self):
        data = [np.array([0.1]), np.array([0.1])]
        result = self.radial_distance_3d.compute(data)
        self.assertIsInstance(result, RadialDistance3D)

    def test_compute_error(self):
        data = [np.array([0.1]), np.array([0.1])]
        errors = [np.array([0.01]), np.array([0.01])]
        result = self.radial_distance_3d.compute_error(data, errors)
        self.assertIsInstance(result, RadialDistance3D)

class TestTotalVelocity(unittest.TestCase):
    def setUp(self):
        self.total_velocity = TotalVelocity()

    def test_compute(self):
        data = [np.array([10.0]), np.array([20.0])]
        result = self.total_velocity.compute(data)
        self.assertIsInstance(result, TotalVelocity)

    def test_compute_error(self):
        data = [np.array([10.0]), np.array([20.0])]
        errors = [np.array([1.0]), np.array([1.0])]
        result = self.total_velocity.compute_error(data, errors)
        self.assertIsInstance(result, TotalVelocity)

class TestEffectivePotential(unittest.TestCase):
    def setUp(self):
        self.effective_potential = EffectivePotential(shell=True)

    def test_compute(self):
        data = [np.array([10.0]), np.array([20.0]), np.array([30.0])]
        result = self.effective_potential.compute(data)
        self.assertIsInstance(result, EffectivePotential)

    def test_compute_error(self):
        data = [np.array([10.0]), np.array([20.0]), np.array([30.0])]
        errors = [np.array([1.0]), np.array([1.0]), np.array([1.0])]
        result = self.effective_potential.compute_error(data, errors)
        self.assertIsInstance(result, EffectivePotential)

if __name__ == '__main__':
    unittest.main()