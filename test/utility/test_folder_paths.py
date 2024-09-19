import unittest
import os

# Save original environment variable
original_gcasdata = os.getenv('GCASDATA')
# Set up the environment variable for testing before importing the module
os.environ['GCASDATA'] = '/mock/path/to/G-GCAS/ggcas'

# Now import the module
from ggcas.utility.folder_paths import CLUSTER_DATA_FOLDER, CLUSTER_MODEL_FOLDER

class TestFolderPaths(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # Save the original environment variable if it exists
        cls.original_gcasdata = original_gcasdata

    @classmethod
    def tearDownClass(cls):
        # Restore the original environment variable
        if cls.original_gcasdata is not None:
            os.environ['GCASDATA'] = cls.original_gcasdata
        else:
            del os.environ['GCASDATA']

    def test_gcas_data_folder(self):
        self.assertEqual(os.getenv('GCASDATA'), '/mock/path/to/G-GCAS/ggcas')

    def test_cluster_data_folder(self):
        expected_path = '/mock/path/to/G-GCAS/ggcas/data/query/CLUSTER1'
        actual_path = CLUSTER_DATA_FOLDER('cluster1')
        print(f"Expected: {expected_path}, Actual: {actual_path}")
        self.assertEqual(actual_path, expected_path)

    def test_cluster_model_folder(self):
        expected_path = '/mock/path/to/G-GCAS/ggcas/data/models/CLUSTER1'
        actual_path = CLUSTER_MODEL_FOLDER('cluster1')
        print(f"Expected: {expected_path}, Actual: {actual_path}")
        self.assertEqual(actual_path, expected_path)

if __name__ == '__main__':
    unittest.main()