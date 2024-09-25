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
from ggcas._utility import folder_paths as fn
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
    tn: str
        Tracking number of the data to search for.
    fold : str
        The name of the folder to search for the tracking number, or the complete
        folder path, including the tn, in case it is a data folder (then 'tn' must
        be none, as default).
    key : str, optional
        A key which identify specific files to return.

    Returns
    -------
    fl : list os str
        List of sorted files inside the folder.

    Raises
    ------
    FileNotFoundError
        If the specified path does not exist.
    TypeError
        If the 'key' argument is not a string.

    Examples
    --------
    Here are some examples regarding the use of the 'key' argument. Let's say
    we need a list of files inside ''tn = '20160516_114916' '' for GC 'ngc104'

        >>> tn = '20160516_114916'
        >>> get_file_list(tn=tn) # if only that [gc_name] folder has that tn folder
        ['.../G-GCAS/ggcas/data/query/[gc_name]/[tn]/query_data.txt',
         '.../G-GCAS/ggcas/data/query/[gc_name]/[tn]/spatial_data.txt',
         '.../G-GCAS/ggcas/data/query/[gc_name]/[tn]/velocity_data.txt',
         '.../G-GCAS/ggcas/data/query/[gc_name]/[tn]/query_info.ini',
         '.../G-GCAS/ggcas/data/query/[gc_name]/[tn]/dynamical_data.txt']

    Let's suppose we want only the list of 'xxx_data.txt' files:

        >>> get_file_list(tn=tn, key='_data')
        ['.../G-GCAS/ggcas/data/query/[gc_name]/[tn]/query_data.txt',
         '.../G-GCAS/ggcas/data/query/[gc_name]/[tn]/spatial_data.txt',
         '.../G-GCAS/ggcas/data/query/[gc_name]/[tn]/velocity_data.txt',
         '.../G-GCAS/ggcas/data/query/[gc_name]/[tn]/dynamical_data.txt']

    This function can be used to retrieve a list of folders too. Say we want to
    know which cluster has king-model data available:

        >>> foldpath = '.../G-GCAS/ggcas/data/models/'
        >>> fold_list = get_file_list(fold=foldpath)
        >>> fold_list
        ['.../G-GCAS/ggcas/data/models/NGC104',
         '.../G-GCAS/ggcas/data/models/NGC4372',
         '.../G-GCAS/ggcas/data/models/NGC6121']

    These are the folders of the clusters which have king-model data available. indeed,
    if we check:
        
        >>> get_file_list(fold=fold_list[0])
        ['.../G-GCAS/ggcas/data/models/NGC104/SM_king.txt']
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

def get_kwargs(names:tuple, default, kwargs):
    """
    Gets a tuple of possible kwargs names for a variable and checks if it was
    passed, and in case returns it.

    Parameters
    ----------
    names : tuple
        Tuple containing all the possible names of a variable which can be passed
        as a **kwargs argument.
    default : any type
        The default value to assign the requested key if it doesn't exist.
    kwargs : dict
        The dictionary of variables passed as 'Other Parameters'.

    Returns
    -------
    key : value of the key
        The value of the searched key if it exists. If not, the default value will
        be returned.
    """
    possible_keys = names
    for key in possible_keys:
        if key in kwargs:
            return kwargs[key]
    return default

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