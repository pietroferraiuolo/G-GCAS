import unittest, mock
from unittest.mock import patch
import os
from ggcas.utility import folder_paths as fn

class TestFolderPaths(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # Save the original environment variable if it exists
        cls.original_gcasdata = os.getenv('GCASDATA')
        # Set up the environment variable for testing
        os.environ['GCASDATA'] = '/mock/path/to/G-GCAS/ggcas'

    @classmethod
    def tearDownClass(cls):
        # Restore the original environment variable
        if cls.original_gcasdata is not None:
            os.environ['GCASDATA'] = cls.original_gcasdata
        else:
            del os.environ['GCASDATA']

    def setUp(self):
        #Reset the environment variable after each test
        os.environ['GCASDATA'] = '/mock/path/to/G-GCAS/ggcas'

    @patch('os.getenv', return_value='/mock/path/to/G-GCAS/ggcas')
    @patch('ggcas.utility.folder_paths.BASE_PATH', return_value='os.getenv("GCASDATA")')
    def test_gcas_data_folder(self, mock_getenv, mock_BASE_PATH):
        self.setUp()
        self.assertEqual(fn.BASE_DATA_PATH, '/mock/path/to/G-GCAS/ggcas/data')
        self.assertEqual(fn.QUERY_DATA_FOLDER, '/mock/path/to/G-GCAS/ggcas/data/query')

    @patch('os.getenv', return_value='/mock/path/to/G-GCAS/ggcas')
    @patch('ggcas.utility.folder_paths.BASE_PATH', return_value='os.getenv("GCASDATA")')
    def test_cluster_data_folder(self, mock_getenv, mock_BASE_PATH):
        expected_path = '/mock/path/to/G-GCAS/ggcas/data/query/CLUSTER1'
        actual_path = fn.CLUSTER_DATA_FOLDER('cluster1')
        print(f"Expected: {expected_path}, Actual: {actual_path}")
        self.assertEqual(actual_path, expected_path)

    @patch('os.getenv', return_value='/mock/path/to/G-GCAS/ggcas')
    @patch('ggcas.utility.folder_paths.BASE_PATH', return_value='os.getenv("GCASDATA")')
    def test_cluster_model_folder(self, mock_getenv, mock_BASE_PATH):
        expected_path = '/mock/path/to/G-GCAS/ggcas/data/models/CLUSTER1'
        actual_path = fn.CLUSTER_MODEL_FOLDER('cluster1')
        print(f"Expected: {expected_path}, Actual: {actual_path}")
        self.assertEqual(actual_path, expected_path)

if __name__ == '__main__':
    unittest.main()