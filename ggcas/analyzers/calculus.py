"""
Author(s):
    - Pietro Ferraiuolo : written in 2024

Copyright : Pietro Ferraiuolo (pietro.ferraiuolo@inaf.it)

Description
-----------

How to Use it
-------------
"""
import os
import subprocess
from typing import Dict, Any
import numpy as np
import sympy as sp
from . import _glpoints
from ggcas.utility import folder_paths as fn, osutils as osu

king_dir = fn.KING_INTEGRATOR_FOLDER
king_exe = os.path.join(king_dir, 'king_integrator')

def gaus_legendre_integrator(fnc, a, b, points):
    """
    Integrates the function fcn(x) between a and b using the Gauss-Legendre method.

    Parameters
    ----------
    fcn : callable
        The function to integrate.
    a : float
        The lower limit of integration.
    b : float
        The upper limit of integration.
    points : int
        The number of points to use for the integration (20, 40, 80, 96).

    Returns
    -------
    area : float
        The integral area of fcn(x) from a to b.

    Notes
    -----
    Note that the imput function must be a 'pythonic' function, and not a sympy
    functions, as it is not implemented yet.
    """
    if points == 20:
        x = np.array(_glpoints.x20)
        w = np.array(_glpoints.w20)
    elif points == 40:
        x = np.array(_glpoints.x40)
        w = np.array(_glpoints.w40)
    elif points == 80:
        x = np.array(_glpoints.x80)
        w = np.array(_glpoints.w80)
    elif points == 96:
        x = np.array(_glpoints.x96)
        w = np.array(_glpoints.w96)
    else:
        raise ValueError("Supported point values are 20, 40, 80, and 96.")
    area = 0.0
    try:
        for i in range(int(points/2)):
            xi = (b - a) / 2.0 * x[i] + (b + a) / 2.0
            area += w[i] * fnc(xi)
            xi = -(b - a) / 2.0 * x[i] + (b + a) / 2.0
            area += w[i] * fnc(xi)
    except TypeError as te:
        raise NotImplementedError(f"{type(fnc)} function type not implemented yet. Define your function using the classic python functions definition") from te
    area *= (b - a) / 2.0
    return area

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
    err_func = error_propagation(func, variables, correlation=corr)
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
        for c in x:
            if c!=0:
                corrs.append(c)
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

def king_integrator(w0, output='profile'):
    r"""
    This function calls a Fortran90 code for the Single-Mass King model integration
    routine.

    Taking as imput a value for the King $W_0$ parameter, it will perform the
    integration of the model, producing a series of output data files, described
    in the 'Returns' section.

    Parameters
    ----------
    w0 : float
        King w0 parameter, which is the central potential well.
    output : str, optional
        Specifies which output file(s) to retain. The default is 'profile'.
        Options :
            - all: All of the below produced files
            - CalCurve: Caloric curve of the system
            - Cv:
            - CvNtK:
            - Er:
            - Etot: total energy of the system
            - params:
            - phi: information about the gravitational potential of the system
            - profiles: (normalized) density and w0 profiles with respect to the dimentionless radial distance from the centre
            - Skin: Surface kinetick energy distribution
            - x0Cv:

    Returns
    -------
    result : str or list
        The full path of the selected output file(s) as string
        (or list of strings if multiple output files have been selected).
    """
    if isinstance(w0, (float, int)):
        w0 = str(w0)
    result = subprocess.run([king_exe, w0], capture_output=True, text=True,
                            cwd=king_dir, check=False)
    if result.returncode != 0:
        print("Error during trhe Fortran90 code execution:")
        print(result.stderr)
    else:
        print("Calling Fortran90 code executor...\n")
        print(result.stdout)
    filelist = osu.get_file_list(fold=king_dir, key='.dat')
    result = []
    if 'all' in output:
        result = filelist
    elif isinstance(output, list):
        for entry in output:
            for i,file in enumerate(filelist):
                if entry in filelist:
                    result.append(file)
                    filelist.pop(i)
    else:
        for i,file in enumerate(filelist):
            if output in file:
                result = file
                filelist.pop(i)
    for file in filelist:
        os.remove(file)
    return result
