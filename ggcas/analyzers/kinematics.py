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
import astropy.units as u

def velocity_conversion(mu, gc_distance, mu_error = 0, gc_distance_error = 0):
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

def density_profile(data):
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