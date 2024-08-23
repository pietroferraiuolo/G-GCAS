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
import numpy as np
import sympy as sp
from astropy import units as u

def angular_separation(ra0, dec0):
    """
    

    Parameters
    ----------
    ra0 : float
        DESCRIPTION.
    dec0 : float
        DESCRIPTION.

    Returns
    -------
    w : sympy expression
        DESCRIPTION.
    variables : list of sympy symbols
        DESCRIPTION.

    """
    ra1, dec1 = sp.symbols('ra1 dec1')
    if isinstance(ra0, u.Quantity):
        ra0 = float(ra0/u.deg)
    if isinstance(dec0, u.Quantity):
        dec0=float(dec0/u.deg)
    d2r = np.pi/180
    w = 2*sp.asin(sp.sqrt(sp.sin((dec0-dec1)*0.5*d2r)**2 + \
                          np.cos(dec0*d2r)*sp.cos(dec1*d2r) * \
                          sp.sin((ra0-ra1)*0.5*d2r)**2))/d2r
    variables = [ra1, dec1]
    func = {'f': w,
            'vars': variables}
    return func

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
