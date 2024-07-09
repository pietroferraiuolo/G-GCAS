"""
Author(s)
---------
    - Pietro Ferraiuolo : Written in 2024

Description
-----------

How to Use
----------

Examples
--------

"""
import os
from astropy.table import QTable
import datetime as dt
datapath    = os.environ['GCASDATA']
querypath   = os.path.join(datapath, 'query')

def load_query(file):
    """
    Loads the data found in the file as an astropy quantity table.

    Parameters
    ----------
    file : str
        Complete file path of the data, obtainable through the ''get_file_list''
        function.

    Returns
    -------
    data : astropy table
        The loaded data of the file.
    """
    data = QTable.read(file, format='ascii.tab')
    return data
        
def get_file_list(gc_name, tn:str=None, key:str=None):
    """
    Returns the file list of a given globular cluster datapath.

    Parameters
    ----------
    gc_name: str
        Name of the globular cluster to search data of.
    key : str, optional
        A key which identify specific files to return.

    Returns
    -------
    fl : list os str
        List of sorted files inside the folder.

    Examples
    --------
    Here are some examples regarding the use of the 'key' argument. Let's say w
    e need a list of files inside ''tn = '20160516_114916' '' for GC 'ngc104' 

        >>> gc_name = 'ngc104'
        >>> tn = '20160516_114916'
        >>> get_file_list(gc_name, tn=tn)
        ['.../G-GCAS/ggcas/data/query/[gc_name]/[tn]/query_data.txt',
         '.../G-GCAS/ggcas/data/query/[gc_name]/[tn]/spatial_data.txt',
         '.../G-GCAS/ggcas/data/query/[gc_name]/[tn]/velocity_data.txt',
         '.../G-GCAS/ggcas/data/query/[gc_name]/[tn]/query_info.ini',
         '.../G-GCAS/ggcas/data/query/[gc_name]/[tn]/dynamical_data.txt']
        
    Let's suppose we want only the list of 'xxx_data.txt' files:
    
        >>> get_file_list(gc_name, tn=tn, key='_data')
        ['.../G-GCAS/ggcas/data/query/[gc_name]/[tn]/query_data.txt',
         '.../G-GCAS/ggcas/data/query/[gc_name]/[tn]/spatial_data.txt',
         '.../G-GCAS/ggcas/data/query/[gc_name]/[tn]/velocity_data.txt',
         '.../G-GCAS/ggcas/data/query/[gc_name]/[tn]/dynamical_data.txt']
    """
    fold = os.path.join(querypath, gc_name)
    if tn is None:
        fl = []
        for item in os.listdir(fold):
            tn = os.path.join(fold, item)
            fl.append(sorted(os.listdir(tn)))
    else:
        fl = sorted(os.listdir(os.pathj.join(fold, tn)))
    if key is not None:
        try:
            selected_list = []
            for file in fl:
                if key in file:
                    selected_list.append(file)
        except TypeError:
            raise TypeError("'key' argument must be a string")
        fl = selected_list
    return fl

def _timestamp():
    """
    Creates a new tracking number in the format 'yyyymmdd_HHMMSS'

    Returns
    -------
    tn : str
        Tracking number.

    """
    tn = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    return tn