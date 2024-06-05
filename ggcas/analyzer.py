"""
Created on May 2024
    -Author: P.Ferraiuolo
"""
import os
import numpy as np
from astropy.table import QTable
import sympy as sp

datapath=os.environ['PYGCASCONF']

def angular_distance(ra0, dec0, ra1, dec1, ra0_err, dec0_err, corr=None):
    '''
    Compute the angular distance between two points on a sphere (coordinates expressed in degrees)
    with it's error, computed through the classic error propagation formula.

    Parameters
    ----------
    ra0 : TYPE
        DESCRIPTION.
    dec0 : TYPE
        DESCRIPTION.
    ra1 : TYPE
        DESCRIPTION.
    dec1 : TYPE
        DESCRIPTION.
    ra0_err : TYPE
        DESCRIPTION.
    dec0_err : TYPE
        DESCRIPTION.
    corr : TYPE
        DESCRIPTION. The default is None

    Returns
    -------
    w : TYPE
        DESCRIPTION.
    we : TYPE
        DESCRIPTION.

    '''
    # Funzioni con sympy
    # Definisci le variabili simboliche
    #ra1, dec1 = sp.symbols('ra1 dec1')
    #ra0, dec0 = sp.symbols('ra0 dec0')
    #d2r = np.pi/180  # Definisci d2r come una variabile simbolica se necessario

    # Usa le funzioni di SymPy per tutte le operazioni matematiche
    #w = 2 * sp.asin(sp.sqrt(sp.sin((dec0 - dec1) * 0.5 * d2r)**2 +
                            #sp.cos(dec0 * d2r) * sp.cos(dec1 * d2r) * sp.sin((ra0 - ra1) * 0.5 * d2r)**2)) / d2r
    
    
    
    
    d2r = np.pi/180  # degrees to radians
    w = 2 * np.arcsin( (np.sin( (dec0-dec1)*0.5 * d2r )**2 +
        np.cos(dec0 * d2r) * np.cos(dec1 * d2r) * np.sin( (ra0-ra1)*0.5 * d2r )**2 )**0.5 ) / d2r

    Omega = np.sin(1/2 * w * d2r)
    denom = 1/np.sqrt(Omega**2 - Omega)

    par_x = denom**2 * np.sin((ra0-ra1)*d2r)*np.cos(dec0*d2r)*np.cos(dec1*d2r)
    par_y = denom**2 * (np.sin((dec0-dec1)*d2r) - 2*np.sin(dec0*d2r)*np.cos(dec1*d2r)*np.sin((ra0-ra1)*0.5*d2r)**2)

    we = np.sqrt((par_x*ra0_err)**2 + (par_y*dec0_err)**2 + par_x*par_y*ra0_err*dec0_err*corr) / d2r

    return w, we

def errPropagation(func, variables, vars_values,  compute=True):
    '''
    

    Parameters
    ----------
    func : sympy expression
        DESCRIPTION.
    variables : list of sympy symbols
        DESCRIPTION.
    errors : float | ArrayLike
        DESCRIPTION.
    corr : float | ArrayLike, optional
        DESCRIPTION. If None, no correlation is assumed and cond is the NxN identity matrix.

    Returns
    -------
    error_propagation : TYPE
        DESCRIPTION.

    '''
    
    corr = np.eye(len(variables))
    errors=[]
    for ii in range(len(variables)):
        errors.append(sp.symbols('epsilon_{}'.format(variables[ii]))) 
        for jj in range(len(variables)):
            if ii != jj:
                corr[ii][jj] = sp.symbols('rho_{}{}'.format(ii, jj))
        
        
    # Calcola le derivate parziali
    partials = [sp.diff(func, var) for var in variables]

    # Primo termine della formula (somma degli errori quadratici)
    sum_of_squares = sum((partials[i]**2 * errors[i]**2) for i in range(len(variables)))

    # Secondo termine della formula (somma delle correlazioni)
    sum_of_correlations = 0
    for ii in range(len(variables)):
        for jj in range(len(variables)):
            if ii != jj:
                sum_of_correlations += (partials[ii] * partials[jj] *
                                        errors[ii] * errors[jj] * corr[ii][jj])

    # Propagazione dell'errore
    error_formula = sp.sqrt(sum_of_squares + 2 * sum_of_correlations)
    
    if compute==True:
        values = dict()
        for xx in range(len(variables)):
            values[variables[xx]] = 0
        
        computed_error = np.zeros(len(vars_values[0]))
        for ii in range(len(vars_values[0])):
            for xx in range(len(variables)):
                values[variables[xx]] = vars_values[xx][ii]
                
            computed_error[ii] = sp.N(error_formula.subs(values))
        
        return computed_error
    else:
        return error_formula

def computeError(err_func, variables, vars_values):
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
    values = dict()
    for x in range(len(variables)):
        values[variables[x]] = 0
    
    computed_error = np.zeros(len(vars_values[0]))
    for i in range(len(vars_values[0])):
        for x in range(len(variables)):
            values[variables[x]] = vars_values[x][i]
            
        computed_error[i] = sp.N(err_func.subs(values))
    
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
    filelist = os.listdir(datapath+name.upper()+'/')
    
    for ii in filelist:
        print(ii)
        
    return filelist