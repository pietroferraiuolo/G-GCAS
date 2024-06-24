"""
Created on May 2024
    -Author: P.Ferraiuolo
"""
from typing import Dict, Any
import numpy as np
import sympy as sp

def errPropagation(func, variables, correlation=False) -> Dict[str, Any]:
    """
    Computes the imput function's error with the standard error propagation
    formula.

    Parameters
    ----------
    func : sympy expression
        DESCRIPTION.
    variables : list of sympy symbols
        DESCRIPTION.
    corr : boolean, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    result : dict
        A dictionary containing the results.
        Keys:
        - "error_formula" : sympy expression
        - "error_variables" : list of sympy symbols
        - "correlations" : list of sympy symbols (only if corr=True)
    """
    if correlation is True:
        corr = np.ones((len(variables), len(variables))) - np.eye(len(variables))
    else:
        corr = np.eye(len(variables))
    errors = []
    for i, var in enumerate(variables):
        errors.append(sp.symbols('epsilon_{}'.format(var)))
        for j, _ in enumerate(variables):
            if i != j:
                corr[i][j] = corr[i][j]*sp.symbols('rho_{}_{}'.format(i, j))
    # Calcola le derivate parziali
    partials = [sp.diff(func, var) for var in variables]
    # Primo termine della formula (somma degli errori quadratici)
    sum_of_squares = sum((partials[i]**2 * errors[i]**2) for i, _ in enumerate(variables))
    # Secondo termine della formula (somma delle correlazioni)
    sum_of_correlations = 0
    for i, _ in enumerate(variables):
        for j, _ in enumerate(variables):
            if i != j:
                sum_of_correlations +=(partials[i]*partials[j]*errors[i]*errors[j]*corr[i][j])
    # Propagazione dell'errore
    error_formula = sp.sqrt(sum_of_squares + 2 * sum_of_correlations)
    returns = {
        "error_formula": error_formula,
        "error_variables": errors
    }
    if correlation:
        returns["correlations"] = corr
    return returns

def angularSeparation():
    """
    

    Returns
    -------
    w : sympy expression
        DESCRIPTION.
    variables : list of sympy symbols
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
    r : sympy expression
        DESCRIPTION.
    parallax : sympy symbol
        DESCRIPTION.

    """
    parallax = sp.symbols('omega')
    r = 1/parallax
    return r, parallax

def radDistance_2D():
    """
    

    Returns
    -------
    r2d : sympy expression
        DESCRIPTION.
    variables : list of sympy symbols
        DESCRIPTION.

    """
    rx, parallax = sp.symbols('r_gc theta')
    r_2d = rx*sp.tan(parallax)
    return r_2d, [rx, parallax]

def radDistance_3D():
    """
    

    Returns
    -------
    R : sympy expression
        DESCRIPTION.
    variables : list of sympy symbols
        DESCRIPTION.

    """
    d, r2d = sp.symbols('d, r_2d')
    r_3d = sp.sqrt(d**2 + r2d**2)
    return r_3d, [d, r2d]
