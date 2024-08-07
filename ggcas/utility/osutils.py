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
from ggcas.utility import folder_paths as fn
datapath    = fn.BASE_DATA_PATH
querypath   = fn.QUERY_DATA_FOLDER

def load_query(file):
    """
    Loads the data found in the file as an astropy quantity table.

    Parameters
    ----------
    file : str
        Complete file path of the data, obtainable through the 'get_file_list'
        function.

    Returns
    -------
    data : astropy table
        The loaded data of the file.
    """
    data = QTable.read(file, format='ascii.tab')
    return data

def get_file_list(tn=None, fold=None, key:str=None):
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
    Here are some examples regarding the use of the 'key' argument. Let's say
    we need a list of files inside ''tn = '20160516_114916' '' for GC 'ngc104'

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
    if tn is None and fold is not None:
        fl = sorted([os.path.join(fold, file) \
                     for file in os.listdir(fold)])
    else:
        if fold is None:
            fold = _findTracknum(tn, complete_path=True)
            fl = sorted([os.path.join(fold, tn, file) \
                         for file in os.listdir(os.path.join(fold, tn))])
        else:
            try:
                paths = _findTracknum(tn, complete_path=True)
                if isinstance(paths, str):
                    paths = [paths]
                for path in paths:
                    if fold in path.split('/')[-2]:
                        fl = sorted([os.path.join(path, file) \
                                                 for file in os.listdir(path)])
                    else:
                        raise Exception
            except Exception as exc:
                raise FileNotFoundError(f"Invalid Path: no data found for '.../{fold}/{tn}'") from exc
    if key is not None:
        try:
            selected_list = []
            for file in fl:
                if key in file.split('/')[-1]:
                    selected_list.append(file)
        except TypeError as err:
            raise TypeError("'key' argument must be a string") from err
        fl = selected_list
    return fl

def tnlist(gc_name:str):
    """


    Parameters
    ----------
    gc_name : str
        DESCRIPTION.

    Returns
    -------
    None.

    """
    basepath = fn.CLUSTER_DATA_FOLDER(gc_name)
    tn = os.listdir(basepath)
    tns = sorted([os.path.join(basepath, tt) for tt in tn])
    return tns

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

def _findTracknum(tn, complete_path:bool=False):

    """
    Search for the tracking number given in input within all the data path subfolders.

    Parameters
    ----------
    tn : str
        Tracking number to be searched.
    complete_path : bool, optional
        Option for wheter to return the list of full paths to the folders which
        contain the tracking number or only their names.

    Returns
    -------
    tn_path : list of str
        List containing all the folders (within the OPTData path) in which the
        tracking number is present, sorted in alphabetical order.

    """
    tn_path = []
    for fold in os.listdir(querypath):
        search_fold = os.path.join(querypath, fold)
        if tn in os.listdir(search_fold):
            if complete_path:
                tn_path.append(os.path.join(search_fold, tn))
            else:
                tn_path.append(fold)
    path_list = sorted(tn_path)
    if len(path_list)==1:
        path_list = path_list[0]
    return path_list

def _get_kwargs(names:tuple, default, **kwargs):
    possible_keys = names
    for key in possible_keys:
        if key in kwargs:
            return kwargs[key]
    return default
