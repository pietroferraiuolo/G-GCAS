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
datapath    = os.environ['PYGCASCONF']
querypath   = os.path.join(datapath, 'query')

def loadQuery(name, tn: str):
    """
    

    Parameters
    ----------
    name : TYPE
        DESCRIPTION.
    tn : str
        DESCRIPTION.

    Raises
    ------
    e
        DESCRIPTION.

    Returns
    -------
    data : TYPE
        DESCRIPTION.

    """
    data = os.path.join(querypath, name)
    file = os.path.join(querypath, (tn+'.txt'))
    try:
        data = QTable.read(file, format='ascii.tab')
        return data
    except FileNotFoundError as e:
        raise e(f"File {file} does not exist.")
        
def dataList(name: str):
    """
    

    Parameters
    ----------
    name : str
        DESCRIPTION.

    Returns
    -------
    filelist : TYPE
        DESCRIPTION.

    """
    filelist = os.listdir(os.path.join(datapath, name.upper()))
    for ii in filelist:
        print(ii)
    return filelist