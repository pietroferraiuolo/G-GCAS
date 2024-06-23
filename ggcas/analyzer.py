"""
Created on May 2024
    -Author: P.Ferraiuolo
"""
import os
import numpy as np
from astropy.table import QTable
import sympy as sp
from ggcas import functions as gfunc

datapath=os.environ['PYGCASCONF']
loaded_name = ''

def computeError(func, variables, variables_values, corr=False):
    '''
    

    Parameters
    ----------
    err_func : TYPE
        DESCRIPTION.
    variables : TYPE
        DESCRIPTION.
    vars_values : TYPE
        DESCRIPTION.

    Returns
    -------
    computed_error : TYPE
        DESCRIPTION.

    '''
    err_func = gfunc.errPropagation(func, variables, correlation=corr)
    values = dict()
    for x in range(len(variables)):
        values[variables[x]] = 0
    computed_error = np.zeros(len(variables_values[0]))
    for i in range(len(variables_values[0])):
        for x in range(len(variables)):
            values[variables[x]] = variables_values[x][i]
        computed_error[i] = sp.N(err_func['error_formula'].subs(values))
    return computed_error
    
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
    global loaded_name
    loaded_name = name    
    filelist = os.listdir(datapath+name.upper()+'/')
    for ii in filelist:
        print(ii)
    return filelist
