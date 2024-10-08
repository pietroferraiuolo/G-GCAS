"""
Module: functions.py

Author(s)
---------
    - Pietro Ferraiuolo : Written in 2024

Description
-----------
This module provides a set of functions to compute various astronomical and
physical quantities such as angular separation, line-of-sight distance, radial
distances in 2D and 3D, total velocity, and effective gravitational potential.

Functions
---------
- angular_separation(ra0=None, dec0=None)
    Computes the angular separation between two points in the sky.

- los_distance()
    Computes the line-of-sight distance based on parallax.

- radial_distance_2d(analytical_w=False, **params)
    Computes the 2D-projected radial distance of a source from the center of a cluster or given RA/DEC coordinates.

- radial_distance_3d(gc_distance=None, analytical_r2d: bool = False, analytical_w: bool = False)
    Computes the 3D radial distance of a source from the center of a cluster or given RA/DEC coordinates.

- total_velocity()
    Computes the total velocity based on the given velocity components.

- effective_potential(shell: bool = False)
    Computes the effective gravitational potential, optionally considering a shell model.

How to Use
----------
1. Import the module:
    ```python
    import functions
    ```

2. Call the desired function with appropriate parameters:
    ```python
    result = functions.angular_separation(ra0=10.0, dec0=20.0)
    ```

Examples
--------
Example usage of `angular_separation` function:
    ```python
    import functions
    from astropy import units as u

    ra0 = 10.0 * u.deg
    dec0 = 20.0 * u.deg
    result = functions.angular_separation(ra0, dec0)
    print(result)
    ```

Example usage of `los_distance` function:
    ```python
    import functions

    result = functions.los_distance()
    print(result)
    ```

Example usage of `radial_distance_2d` function:
    ```python
    import functions

    result = functions.radial_distance_2d(analytical_w=True, ra0=10.0, dec0=20.0)
    print(result)
    ```

Example usage of `radial_distance_3d` function:
    ```python
    import functions
    from astropy import units as u

    gc_distance = 1000 * u.pc
    result = functions.radial_distance_3d(gc_distance=gc_distance, analytical_r2d=True, analytical_w=True)
    print(result)
    ```

Example usage of `total_velocity` function:
    ```python
    import functions

    result = functions.total_velocity()
    print(result)
    ```

Example usage of `effective_potential` function:
    ```python
    import functions

    result = functions.effective_potential(shell=True)
    print(result)
    ```
"""
import numpy as np
import sympy as sp
from astropy import units as u

def angular_separation(ra0=None, dec0=None):
    """


    Parameters
    ----------
    ra0 : float, optional
        DESCRIPTION.
    dec0 : float, optional
        DESCRIPTION.

    Returns
    -------
    func : dict
        DESCRIPTION.

    """
    ra1, dec1 = sp.symbols("alpha_1 \delta_1")
    variables = [ra1, dec1]
    d2r = np.pi / 180
    if ra0 is not None and dec0 is not None:
        if isinstance(ra0, u.Quantity):
            ra0 = float(ra0 / u.deg)
        if isinstance(dec0, u.Quantity):
            dec0 = float(dec0 / u.deg)
        costerm = np.cos(dec0 * d2r)
    else:
        ra0, dec0 = sp.symbols("alpha_0 \delta_0")
        variables.append(ra0)
        variables.append(dec0)
        costerm = sp.cos(dec0 * d2r)
    w = (
        2
        * sp.asin(
            sp.sqrt(
                sp.sin((dec0 - dec1) * 0.5 * d2r) ** 2
                + costerm * sp.cos(dec1 * d2r) * sp.sin((ra0 - ra1) * 0.5 * d2r) ** 2
            )
        )
        / d2r
    )
    func = {"f": w, "vars": variables}
    return func

def los_distance():
    """


    Returns
    -------
    func : TYPE
        DESCRIPTION.

    """
    parallax = sp.symbols("omega")
    r = 1 / parallax
    func = {"f": r, "vars": parallax}
    return func

def radial_distance_2d(analytical_w=False, **params):
    """
    Returns the formula for the 2d-projection on the plane of the sky of the radial
    distance of a source from the center of the cluster (or from the given RA/DEC
    coordinates).

    Parameters
    ----------
    analytical_w : bool, optional
        Whether to have the full expression of the angular separation in the returned
        formula, or the representing symbol (if, e.g, has been altrady calculated).
        The default is False.

    Other Parameters
    ----------------
    **params : dict
        Additional parameters, callback for the 'angular_separation' function.
        ra0 : float
            RA coordinate from which compute the angular separation.
        dec0 : float
            DEC coordinate from which compute the angular separation.

    Returns
    -------
    func : dict
        Dictionary containing the formula and the relative variables for the
        2D-projected radial distance.

    """
    rgc = sp.symbols("r_gc")
    variables = [rgc]
    if analytical_w is False:
        w = sp.symbols("theta_x0")
        variables.append(w)
    else:
        ww = angular_separation()
        for var in ww["vars"]:
            variables.append(var)
        w = ww["f"]
    r_2d = rgc * sp.tan(w)
    func = {"f": r_2d, "vars": variables}
    return func

def radial_distance_3d(
    gc_distance=None, analytical_r2d: bool = False, analytical_w: bool = False
):
    """


    Returns
    -------
    func : dict
        Dictionary containing the formula and the relative variables for the
        3D radial distance.

    Notes
    -----
    If gc_distance is given as argument, the computation of this formula returns
    the radial distance in parsec, otherwise the unit will be the one used when
    passing the cluster's distance manually in the computation.
    """
    rx = los_distance()
    variables = [rx["vars"]]
    if gc_distance is not None:
        rgc = (
            gc_distance.to(u.pc).value
            if isinstance(gc_distance, u.Quantity)
            else gc_distance
        )
    else:
        rgc = sp.symbols("r_gc")
        variables.append(rgc)
    if (analytical_r2d, analytical_w) == (False, True):
        print("WARNING: 'analytical_r2d is False, no effect for analytical_w")
    if analytical_r2d:
        r2_d = radial_distance_2d(analytical_w=analytical_w)
        r2d = r2_d["f"]
        for v in r2_d["vars"]:
            variables.append(v)
        if gc_distance is not None:
            r2d = sp.N(r2d.subs({rgc: gc_distance}))
        variables.remove(sp.symbols("r_gc"))
    else:
        r2d = sp.symbols("r_2d")
        variables.append(r2d)
    D = rx["f"] - rgc
    r_3d = sp.sqrt(D**2 + r2d**2)
    func = {"f": r_3d, "vars": variables, "d": D}
    return func

def total_velocity():
    """


    Returns
    -------
    func : dict
        Dictionary containing the formula and the relative variables for the
        total velocity.
    """
    vx, vy = sp.symbols("v_x, v_y")
    V = 1.5 * (vx**2 + vy**2)
    variables = [vx, vy]
    func = {"f": V, "vars": variables}
    return func

def effective_potential(shell: bool = False):
    """


    Parameters
    ----------
    shell : bool, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    func : dict
        Dictionary containing the formula and the relative variables for the
        effective gravitational potential.
    """
    if shell:
        dN, x, lnB = sp.symbols("\Delta\ N, x, lnB")
        A, a, b, r, m, sigma, dr, dx = sp.symbols(
            "A, alpha, beta, r_*, m, sigma, dr, dx"
        )
        lnb = (
            16 * sp.sqrt(2) * A * (sp.pi * b * r) ** 2 * a * (m * sigma) ** 3 * dr * dx
        )
        variables = {"x": [dN, x, lnB], "const": [lnb, A, a, b, r, m, sigma, dr, dx]}
        poteff = lnB - sp.ln(dN / sp.sqrt(x)) - x
    else:
        x, w = sp.symbols("x, w")
        variables = [x, w]
        poteff = -sp.ln(1 - sp.exp(x - w))
    func = {"f": poteff, "vars": variables}
    return func
