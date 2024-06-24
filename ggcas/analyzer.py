"""
Created on May 2024
    -Author: P.Ferraiuolo
"""
import os
import numpy as np
import sympy as sp
import astropy.units as u
from astropy.table import QTable
from ggcas import functions as gfunc

datapath    = os.environ['PYGCASCONF']
querypath   = os.path.join(datapath, 'query')

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

    '''#!!! to check
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

def velocityConversion(mu, gc_distance, mu_error = 0, gc_distance_error = 0):
    """
    Converts the proper motion into velocities in km/s, with its error, if provided.

    Parameters
    ----------
    mu : TYPE
        DESCRIPTION.
    gc_distance : TYPE
        DESCRIPTION.
    mu_error : TYPE, optional
        DESCRIPTION. The default is 0.
    gc_distance_error : TYPE, optional
        DESCRIPTION. The default is 0.

    Returns
    -------
    vkms : TYPE
        DESCRIPTION.
    vkms_err : TYPE
        DESCRIPTION.

    """
    vkms = mu.to(u.mas/u.yr).to(u.rad/u.s)*gc_distance.to(u.kpc).to(u.km) / u.rad
    vkms_err = np.sqrt(gc_distance**2 * mu_error.to(u.rad/u.s)**2 + \
                       mu.to(u.rad/u.s)**2 * gc_distance_error**2)/u.rad
    return vkms, vkms_err

def densityProfile(data):
    """
    

    Parameters
    ----------
    data : TYPE
        DESCRIPTION.

    Returns
    -------
    rho : TYPE
        DESCRIPTION.

    """
    n_bin = int(1.5*len(data)**0.5)
    dh = np.histogram(data, bins=n_bin)
    V = np.zeros(n_bin)
    rr1 = np.zeros(n_bin+1)
    rr2 = np.zeros(n_bin+1)
    bw = data.max()/n_bin
    rr1[0] = 0.
    rr2[0] = bw
    for x in range (0, n_bin):
      V[x] = (4/3)*np.pi*(rr2[x]**3 - rr1[x]**3)
      rr1[x+1] = rr1[x] + bw
      rr2[x+1] = rr2[x] + bw
    rho = dh[0]/(V**3)
    return rho
    
def loadQuery(name, tn: str):
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
    data = os.path.join(querypath, name)
    file = os.path.join(querypath, (tn+'.txt'))
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
    filelist = os.listdir(os.path.join(datapath, name.upper()))
    for ii in filelist:
        print(ii)
    return filelist