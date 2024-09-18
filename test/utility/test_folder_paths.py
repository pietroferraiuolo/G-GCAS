import unittest
import os
from ggcas.utility.folder_paths import CLUSTER_DATA_FOLDER, CLUSTER_MODEL_FOLDER

class TestFolderPaths(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # Set up the environment variable for testing
        os.environ['GCASDATA'] = '/mock/path/to/G-GCAS/ggcas'

    def test_cluster_data_folder(self):
        expected_path = '/mock/path/to/G-GCAS/ggcas/data/query/CLUSTER1'
        self.assertEqual(CLUSTER_DATA_FOLDER('cluster1'), expected_path)

    def test_cluster_model_folder(self):
        expected_path = '/mock/path/to/G-GCAS/ggcas/data/models/CLUSTER1'
        self.assertEqual(CLUSTER_MODEL_FOLDER('cluster1'), expected_path)

if __name__ == '__main__':
    unittest.main()