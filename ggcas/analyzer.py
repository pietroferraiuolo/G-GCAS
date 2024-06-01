"""
Created on May 2024
    -Author: P.Ferraiuolo
"""
import os
import numpy as np
from astropy.table import QTable

datapath=os.environ['PYGCASCONF']

def loadQuery(tn: str):
    '''
    

    Parameters
    ----------
    tn : str
        DESCRIPTION.

    Returns
    -------
    data : TYPE
        DESCRIPTION.

    '''
    file = os.path.join(datapath, (tn+'.txt'))
    
    try:
        data = QTable.read(file, format='ascii.tab')
        return data
    except FileNotFoundError as e:
        raise e
        
def dataList(name: str):
    '''
    

    Parameters
    ----------
    name : str
        DESCRIPTION.

    Returns
    -------
    filelist : TYPE
        DESCRIPTION.

    '''
    filelist = os.listdir(datapath+name.upper()+'/')
    
    for ii in filelist:
        print(ii)
        
    return filelist
