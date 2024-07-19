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
from typing import Dict, Any
import numpy as np
import sympy as sp

def error_propagation(func, variables, correlation:bool=False) -> Dict[str, Any]:
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
    if correlation:
        corr = np.ones((len(variables), len(variables)), dtype=object) -\
                                           np.eye(len(variables), dtype=object)
    else:
        corr = np.eye(len(variables), dtype=object)
    errors = []
    for i, var1 in enumerate(variables):
        errors.append(sp.symbols('epsilon_{}'.format(var1)))
        for j, var2 in enumerate(variables):
            if i != j:
                corr[i][j] = corr[i][j]*sp.symbols(f'rho_{var1}_{var2}')
    # Partial derivatives computation
    partials = [sp.diff(func, var) for var in variables]
    # Quadratic errors sum
    sum_of_squares = sum((partials[i]**2 * errors[i]**2) for i, _ in enumerate(variables))
    # Correlation sums
    sum_of_correlations = 0
    for i, _ in enumerate(variables):
        for j, _ in enumerate(variables):
            if i != j:
                sum_of_correlations +=(partials[i]*partials[j]*errors[i]*errors[j]*corr[i][j])
    # Total error propagation
    error_formula = sp.sqrt(sum_of_squares + sum_of_correlations)
    corrs = []
    for x in corr:
        if x!=0:
            corrs.append(corr)
    returns = {
        "error_formula": error_formula,
        "error_variables": {
            'variables': variables,
            'errors': errors
            }
    }
    if correlation:
        returns["correlations"] = corr
    return returns

def angular_separation():
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
                          sp.cos(dec0*d2r)*sp.cos(dec1*d2r) * \
                          sp.sin((ra0-ra1)*0.5*d2r)**2))/d2r
    variables = [ra1, dec1, ra0, dec0]
    return w, variables

def los_distance():
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

def radial_distance_2d():
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
    variables = [rx, parallax]
    return r_2d, variables


def radial_distance_3d():
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
    variables = [d, r2d]
    return r_3d, variables

def total_velocity():
    """
    

    Returns
    -------
    V : TYPE
        DESCRIPTION.
    list
        DESCRIPTION.

    """
    vx, vy = sp.symbols('v_x, v_y')
    V = 1.5*(vx**2 + vy**2)
    variables = [vx,vy]
    return V, variables

def effective_potential(shell = False):
    """
    

    Parameters
    ----------
    shell : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    poteff : TYPE
        DESCRIPTION.
    variables : TYPE
        DESCRIPTION.

    """
    if shell:
        dN, x = sp.symbols('DeltaN, x')
        A, a, b, r, m, sigma, dr, dx = sp.symbols('A, alpha, beta, r_*, m, sigma, dr, dx')
        lnB = 16*sp.sqrt(2)*A*(sp.pi*b*r)**2 * a*(m*sigma)**3 * dr*dx
        constants = [lnB, [A, a, b, r, m, sigma, dr, dx]]
        variables = [dN, x]
        variables.append(constants)
        poteff = lnB - sp.ln(dN/sp.sqrt(x)) - x
    else:
        x, w = sp.symbols('x, w')
        variables = [x,w]
        poteff = -sp.ln(1-sp.exp(x-w))
    return poteff, variables
