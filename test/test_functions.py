import unittest
import numpy as np
from astropy import units as u
from grasp import functions as f

class TestAngularSeparation(unittest.TestCase):
    def setUp(self):
        self.ra0 = 10.0 * u.deg
        self.dec0 = 20.0 * u.deg
        self.angular_separation = f.AngularSeparation(self.ra0, self.dec0)

    def test_compute(self):
        data = [np.array([15.0]), np.array([25.0])]
        result = self.angular_separation.compute(data)
        self.assertIsInstance(result, f.AngularSeparation)
        self.assertNotEqual(result.computed_values, None)
        #np.testing.assert_array_almost_equal(result.values, np.array([6.80561008])) 

    def test_compute_error(self):
        data = [np.array([15.0]), np.array([25.0])]
        errors = [np.array([0.1]), np.array([0.1])]
        corr = [np.array([0.1])]
        result = self.angular_separation.compute(data, errors, corr)
        self.assertIsInstance(result, f.AngularSeparation)
        self.assertNotEqual(result.computed_errors, None)
        #np.testing.assert_array_almost_equal(result.errors, np.array([0.09563928]))

class TestLosDistance(unittest.TestCase):
    def setUp(self):
        self.los_distance = f.LosDistance()

    def test_compute(self):
        data = [np.array([0.1])]
        result = self.los_distance.compute(data)
        self.assertIsInstance(result, f.LosDistance)
        self.assertNotEqual(result.computed_values, None)
        #np.testing.assert_array_almost_equal(result.values, np.array([10.]))

    def test_compute_error(self):
        data = [np.array([0.1])]
        errors = [np.array([0.01])]
        result = self.los_distance.compute(data, errors)
        self.assertIsInstance(result, f.LosDistance)
        self.assertNotEqual(result.computed_errors, None)
        #np.testing.assert_array_almost_equal(result.errors, np.array([1.]))

class TestRadialDistance2D(unittest.TestCase):
    def setUp(self):
        self.gc_distance = 1000 * u.pc
        self.radial_distance_2d = f.RadialDistance2D(self.gc_distance)

    def test_compute(self):
        data = [np.array([0.1])]
        result = self.radial_distance_2d.compute(data)
        self.assertIsInstance(result, f.RadialDistance2D)
        self.assertNotEqual(result.computed_values, None)

    def test_compute_error(self):
        data = [np.array([1., 2])]
        errors = [np.array([0.01, 0.02])]
        result = self.radial_distance_2d.compute(data, errors)
        self.assertIsInstance(result, f.RadialDistance2D)
        assert result.computed_errors is not None
        #self.assertNotEqual(result.computed_errors, None)

class TestRadialDistance3D(unittest.TestCase):
    def setUp(self):
        self.gc_distance = 1851
        self.r3d = f.RadialDistance3D()
        self.r3d_ext = f.RadialDistance3D(gc_distance=self.gc_distance)

    def test_compute(self):
        data = [np.array([0.1]), np.array([0.1])]
        result = self.r3d.compute(data)
        self.assertIsInstance(result, f.RadialDistance3D)
        self.assertNotEqual(result.computed_values, None)

    def test_compute_error(self):
        data = [np.array([0.1]), np.array([0.1])]
        errors = [np.array([0.01]), np.array([0.01])]
        result = self.r3d_ext.compute(data, errors)
        self.assertIsInstance(result, f.RadialDistance3D)
        self.assertNotEqual(result.computed_errors, None)

class TestTotalVelocity(unittest.TestCase):
    def setUp(self):
        self.total_velocity = f.TotalVelocity()

    def test_compute(self):
        data = [np.array([10.0]), np.array([20.0])]
        result = self.total_velocity.compute(data)
        self.assertIsInstance(result, f.TotalVelocity)
        self.assertNotEqual(result.computed_values, None)

    def test_compute_error(self):
        data = [np.array([10.0]), np.array([20.0])]
        errors = [np.array([1.0]), np.array([1.0])]
        result = self.total_velocity.compute(data, errors)
        self.assertIsInstance(result, f.TotalVelocity)
        self.assertNotEqual(result.computed_errors, None)

class TestEffectivePotential(unittest.TestCase):
    def setUp(self):
        self.effective_potential = f.EffectivePotential(shell=True)

    def test_compute(self):
        data = [np.array([10.0]), np.array([20.0]), np.array([30.0])]
        result = self.effective_potential.compute(data)
        self.assertIsInstance(result, f.EffectivePotential)
        self.assertNotEqual(result.computed_values, None)

    def test_compute_error(self):
        data = [np.array([10.0]), np.array([20.0]), np.array([30.0])]
        errors = [np.array([1.0]), np.array([1.0]), np.array([1.0])]
        result = self.effective_potential.compute(data, errors)
        self.assertIsInstance(result, f.EffectivePotential)
        self.assertNotEqual(result.computed_errors, None)

if __name__ == '__main__':
    unittest.main()