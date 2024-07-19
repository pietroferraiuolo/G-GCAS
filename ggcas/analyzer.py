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
import astropy.units as u
from ggcas import functions as gfunc

def compute_error(func, variables, var_data, var_errors, corr:bool=False,
                                                        corr_values:list=None):
    """


    Parameters
    ----------
    func : TYPE
        DESCRIPTION.
    variables : TYPE
        DESCRIPTION.
    var_data : TYPE
        DESCRIPTION.
    var_errors : TYPE
        DESCRIPTION.
    corr : bool, optional
        DESCRIPTION. The default is False.
    corr_values : list, optional
        DESCRIPTION. The default is None.

    Returns
    -------
    computed_error : TYPE
        DESCRIPTION.

    """
    err_func = gfunc.error_propagation(func, variables, correlation=corr)
    func = err_func['error_formula']
    errors = err_func['error_variables']['errors']
    vars_to_pass = []
    vals_to_pass = []
    for i, var in enumerate(variables):
        vars_to_pass.append(var)
        vals_to_pass.append(var_data[i])
    for i,_  in enumerate(variables):
        vars_to_pass.append(errors[i])
        vals_to_pass.append(var_errors[i])
    if corr:
        corr_values = [corr_values]
        N = len(corr_values)
        if N > 1:
            for i in range(N):
                vars_to_pass.append(err_func['correlations'][i])
                vals_to_pass.append(corr_values[i])
        else:
            vars_to_pass.append(err_func['correlations'])
            vals_to_pass.append(corr_values[0])
    computed_error = compute_numerical_function(func, vars_to_pass, vals_to_pass)
    return computed_error

def compute_numerical_function(func, variables, var_data):
    """
    Compute the numerical value of a function, passing by the function, it's
    variables and the data associated to the variables.

    Parameters
    ----------
    func : sympy.core.function
        Function to compute. Must be a sympy expression.
    variables : list of sympy variables
        Variables of the function, as sympy symbols, organizaed in a list.
    var_data : list of ndarray
        Numerical values of the variables, organized in a list ordered the same
        way as the order of the variable list.

    Returns
    -------
    computed_func : float | ArrayLike
        List of values of the function computed for each data point.

    """
    val_dicts = []
    for n in range(len(var_data[0])):
        data_dict = {}
        for i, var in enumerate(variables):
            var_name = f"{var}"
            data_dict[var_name] = var_data[i][n]
        val_dicts.append(data_dict)
    computed_func = [float(sp.N(func.subs(vals))) for vals in val_dicts]
    return computed_func

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
