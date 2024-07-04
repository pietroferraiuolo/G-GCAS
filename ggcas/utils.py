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
datapath    = os.environ['PYGCASCONF']
querypath   = os.path.join(datapath, 'query')

def findQuery(tn):
    """
    

    Parameters
    ----------
    tn : TYPE
        DESCRIPTION.

    Returns
    -------
    data : TYPE
        DESCRIPTION.

    """
    data = []
    for fold in os.listdir(querypath):
        cluster_fold = os.path.join(querypath, fold)
        if tn in os.listdir(cluster_fold):
            data.append([cluster_fold, tn])
    return data

def loadQuery(tn: str, cluster_name:str=None, file_extension:str='txt'):
    """
    Loads a query file into an astropy.table QTable (quantity table).

    Parameters
    ----------
    tn : str
        DESCRIPTION.
    cluster_name : str, optional
        DESCRIPTION. The default is None.
    file_extension : str, optional
        DESCRIPTION. The default is 'txt'.

    Raises
    ------
    FileNotFoundError
        DESCRIPTION.

    Returns
    -------
    data : TYPE
        DESCRIPTION.

    """
    ext = '.'+file_extension
    name = tn+ext
    if cluster_name is not None:
        cluster_name = cluster_name.upper()
        data = os.path.join(querypath, cluster_name)
        file = os.path.join(querypath, name)
        try:
            data = QTable.read(file, format='ascii.tab')
        except FileNotFoundError as e:
            raise FileNotFoundError(f"No such file: '{file}'") from e
    else:
        clist = os.listdir(querypath)
        file_found = False
        if len(clist) != 0:
            for fold in clist:
                cluster_fold = os.path.join(querypath, fold)
                if name in os.listdir(cluster_fold):
                    file = os.path.join(cluster_fold, name)
                    data = QTable.read(file, format='ascii.tab')
                    file_found = f'File found in {cluster_fold}'
                    print(file_found)
            if file_found is False:
                raise FileNotFoundError(f"No such file found: '{name}'")
        else:
            raise FileNotFoundError(f"No data found in '{querypath}'")
    return data
        
def dataList(name: str):
    """
    Function that returns the list of files inside a cluster's data folder.

    Parameters
    ----------
    name : str
        Name of the cluster. 
        Each cluster with save data will have an upper case folder with its nam
        e in it. simply pass in the name of the cluster.

    Returns
    -------
    filelist : list of str
        List of files inside the cluster's folder.

    """
    filelist = sorted(os.listdir(os.path.join(datapath, name.upper())))
    for ii in filelist:
        print(ii)
    return filelist

def timestamp():
    """
    Creates a new tracking number in the format 'yyyymmdd_HHMMSS'

    Returns
    -------
    tn : str
        Tracking number.

    """
    tn = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    return tn