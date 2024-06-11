"""
Created on May 2024
    -Author: P.Ferraiuolo
"""
import numpy as np
import sympy as sp

def errPropagation(func, variables, corr=False):
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
    if corr==True:
        corr = np.unos((len(variables),len(variables)))-np.eye(len(variables))
    else:
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

    return error_formula

def angularDistance():
    ra1, ra0, dec1, dec0 = sp.symbols('ra1 ra0 dec1 dec0')
    d2r = np.pi/180
    w = 2*sp.asin(sp.sqrt(sp.sin((dec0-dec1)*0.5*d2r)**2 + sp.cos(dec0*d2r)*sp.cos(dec1*d2r)*sp.sin((ra0-ra1)*0.5*d2r)**2))/d2r
    return w, [ra1, dec1, ra0, dec0]

def losDistance():
    ww = sp.symbols('omega')
    r = 1/ww
    return r, ww

def radDistance_2D():
    rx, ww = sp.symbols('r_gc theta')
    r2d = rx*sp.tan(ww)
    return r2d, [rx, ww]

def radDistance_3D():
    d, r2d = sp.symbols('d, r_2d')
    R = sp.sqrt(d**2 + r2d**2)
    return R, [d, r2d]
