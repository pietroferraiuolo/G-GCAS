import unittest
import os
from unittest.mock import patch
from astropy.table import QTable
from ggcas.utility import osutils as osu
from ggcas.utility import folder_paths as fn
import numpy as np

class TestOsutils(unittest.TestCase):

    @patch('astropy.table.QTable.read')
    def test_load_query(self, mock_read):
        # Mock the QTable.read function
        mock_table = QTable(names=('col1', 'col2'), rows=[(1, 2), (3, 4)])
        mock_read.return_value = mock_table
        # Test the load_query function
        result = osu.load_query('dummy_file')
        self.assertTrue(np.array_equal(result['col1'], mock_table['col1']))
        self.assertTrue(np.array_equal(result['col2'], mock_table['col2']))
        mock_read.assert_called_once_with('dummy_file', format='ascii.tab')

    def test_get_kwargs(self):
        # Test with key present in kwargs
        kwargs = {'key1': 'value1', 'key2': 'value2'}
        result = osu.get_kwargs(('key1', 'key3'), 'default', kwargs)
        self.assertEqual(result, 'value1')
        # Test with key not present in kwargs
        result = osu.get_kwargs(('key3', 'key4'), 'default', kwargs)
        self.assertEqual(result, 'default')

    @patch('os.listdir')
    @patch('os.path.join')
    @patch('ggcas.utility.osutils._findTracknum')
    def test_get_file_list(self, mock_findTracknum, mock_join, mock_listdir):
        # Mock the listdir and join functions
        mock_listdir.return_value = ['file1.txt', 'file2_data.txt', 'file3.txt']
        mock_join.side_effect = lambda *args: '/'.join(args)
        mock_findTracknum.return_value = '/basepath/folder1'

        # Test with tn and key
        result = osu.get_file_list(tn='tn1', key='_data')
        expected = ['/basepath/folder1/tn1/file2_data.txt']
        self.assertEqual(result, expected)

        # Test with fold and no key
        result = osu.get_file_list(fold='/basepath/folder1')
        expected = ['/basepath/folder1/file1.txt', '/basepath/folder1/file2_data.txt', '/basepath/folder1/file3.txt']
        self.assertEqual(result, expected)

    @patch('os.listdir')
    @patch('ggcas.utility.osutils.fn.CLUSTER_DATA_FOLDER')
    def test_tnlist(self, mock_cluster_data_folder, mock_listdir):
        # Mock the CLUSTER_DATA_FOLDER and listdir functions
        mock_cluster_data_folder.return_value = '/basepath/gc_name'
        mock_listdir.return_value = ['tn1', 'tn2', 'tn3']

        # Test the tnlist function
        result = osu.tnlist('gc_name')
        expected = ['/basepath/gc_name/tn1', '/basepath/gc_name/tn2', '/basepath/gc_name/tn3']
        self.assertEqual(result, expected)

    def test_timestamp(self):
        # Test the _timestamp function
        result = osu._timestamp()
        # Check if the result matches the expected format
        self.assertRegex(result, r'\d{8}_\d{6}')

    @patch('os.listdir')
    @patch('os.path.join')
    def test_findTracknum(self, mock_join, mock_listdir):
        # Mock the listdir and join functions
        mock_listdir.return_value = ['folder1', 'folder2', 'folder3']
        mock_join.side_effect = lambda basepath, folder: f"{basepath}/{folder}"

        # Test with complete_path=False
        result = osu._findTracknum('folder1', complete_path=False)
        expected = ['folder1']
        self.assertEqual(result, expected)

        # Test with complete_path=True
        result = osu._findTracknum('folder1', complete_path=True)
        expected = ['/basepath/folder1']
        self.assertEqual(result, expected)

if __name__ == '__main__':
    unittest.main()