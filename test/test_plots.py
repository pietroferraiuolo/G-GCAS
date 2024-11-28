import unittest
import numpy as np
import pandas as pd
from ggcas import plots as gplt

class TestPlots(unittest.TestCase):

    def setUp(self):
        # Sample data for testing
        self.x = np.random.randn(100)
        self.y = np.random.randn(100)
        self.g = np.random.uniform(10, 20, 100)
        self.b_r = np.random.uniform(0, 2, 100)
        self.teff_gspphot = np.random.uniform(3000, 8000, 100)
        self.pmra = np.random.uniform(-5, 5, 100)
        self.pmdec = np.random.uniform(-5, 5, 100)
        self.ra = np.random.uniform(0, 360, 100)
        self.dec = np.random.uniform(-90, 90, 100)
        self.data = np.random.randn(100)
        self.dataerr = np.random.uniform(0.1, 0.5, 100)
        self.sample = pd.DataFrame({'ra': self.ra, 'dec': self.dec, 'pmra': self.pmra, 'pmdec': self.pmdec})

    def test_scatter_2hist(self):
        try:
            gplt.scatter_2hist(self.x, self.y, kde=True, xlabel='X Axis', ylabel='Y Axis', title='Scatter 2Hist')
        except Exception as e:
            self.fail(f"scatter_2hist raised an exception: {e}")

    def test_colorMagnitude(self):
        try:
            gplt.colorMagnitude(self.g, self.b_r, self.teff_gspphot, alpha=0.5, cmap='viridis')
        except Exception as e:
            self.fail(f"colorMagnitude raised an exception: {e}")

    def test_properMotion(self):
        try:
            gplt.properMotion(self.sample, color='blue', s=10, alpha=0.7)
        except Exception as e:
            self.fail(f"properMotion raised an exception: {e}")

    def test_spatial(self):
        try:
            gplt.spatial(self.sample, color='red', s=15, alpha=0.6)
        except Exception as e:
            self.fail(f"spatial raised an exception: {e}")

    def test_histogram(self):
        try:
            result = gplt.histogram(self.data, kde=True, xlabel='Data', out=True)
            self.assertIn('h', result)
            self.assertIn('kde', result)
        except Exception as e:
            self.fail(f"histogram raised an exception: {e}")

    def test_scat_xhist(self):
        try:
            result = gplt.scat_xhist(self.x, self.y, xerr=0.1, xlabel='X Axis', ylabel='Y Axis')
            self.assertEqual(len(result), 2)
        except Exception as e:
            self.fail(f"scat_xhist raised an exception: {e}")

    def test_errorbar(self):
        try:
            gplt.errorbar(self.data, self.dataerr, x=self.x, xerr=0.1, fmt='o', color='green', ecolor='black')
        except Exception as e:
            self.fail(f"errorbar raised an exception: {e}")

if __name__ == '__main__':
    unittest.main()