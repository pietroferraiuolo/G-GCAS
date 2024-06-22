"""
Created on May 2024
    -Author: P.Ferraiuolo
"""
import numpy as np
import sympy as sp

def errPropagation(func, variables, corr=False):
    """
    
    
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
    error_propagation : sympy expression
        DESCRIPTION.
    errors : list of sympy symbols
        DESCRIPTION.
    corr : list of sympy symbols
        DESCRIPTION.
    """    
    if corr is True:
        corr = np.ones((len(variables), len(variables))) - np.eye(len(variables))
    else:
        corr = np.eye(len(variables))
    errors = []
    for i, var in enumerate(variables):
        errors.append(sp.symbols('epsilon_{}'.format(var)))
        for j, _ in enumerate(variables):
            if i != j:
                corr[i][j] = sp.symbols('rho_{}_{}'.format(i, j))
    # Calcola le derivate parziali
    partials = [sp.diff(func, var) for var in variables]
    # Primo termine della formula (somma degli errori quadratici)
    sum_of_squares = sum((partials[i]**2 * errors[i]**2) for i, _ in enumerate(variables))
    # Secondo termine della formula (somma delle correlazioni)
    sum_of_correlations = 0
    for i, _ in enumerate(variables):
        for j, _ in enumerate(variables):
            if i != j:
                sum_of_correlations += (partials[i] * partials[j] * errors[i] * errors[j] * corr[i][j])
    # Propagazione dell'errore
    error_formula = sp.sqrt(sum_of_squares + 2 * sum_of_correlations)
    return error_formula, errors, corr

def angularDistance():
    """
    

    Returns
    -------
    w : TYPE
        DESCRIPTION.
    list
        DESCRIPTION.

    """
    ra1, ra0, dec1, dec0 = sp.symbols('ra1 ra0 dec1 dec0')
    d2r = np.pi/180
    w = 2*sp.asin(sp.sqrt(sp.sin((dec0-dec1)*0.5*d2r)**2 + \
                          sp.cos(dec0*d2r)*sp.cos(dec1*d2r)*sp.sin((ra0-ra1)*0.5*d2r)**2))/d2r
    return w, [ra1, dec1, ra0, dec0]

def losDistance():
    """
    

    Returns
    -------
    r : TYPE
        DESCRIPTION.
    ww : TYPE
        DESCRIPTION.

    """
    ww = sp.symbols('omega')
    r = 1/ww
    return r, ww

def radDistance_2D():
    """
    

    Returns
    -------
    r2d : TYPE
        DESCRIPTION.
    list
        DESCRIPTION.

    """
    rx, ww = sp.symbols('r_gc theta')
    r_2d = rx*sp.tan(ww)
    return r_2d, [rx, ww]

def radDistance_3D():
    """
    

    Returns
    -------
    R : TYPE
        DESCRIPTION.
    list
        DESCRIPTION.

    """
    d, r2d = sp.symbols('d, r_2d')
    r_3d = sp.sqrt(d**2 + r2d**2)
    return r_3d, [d, r2d]
