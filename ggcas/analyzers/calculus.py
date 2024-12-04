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
import multiprocessing as mp
from typing import Dict, Any, List
import numpy as np
from numpy.typing import ArrayLike
import sympy as sp
from ggcas._utility import KING_INTEGRATOR_FOLDER, get_file_list
from . import _glpoints

_king_dir = KING_INTEGRATOR_FOLDER
_king_exe = os.path.join(_king_dir, 'king_integrator')

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

    Notes
    -----
    When the computation exceeds 30 seconds running on all cores, which should
    happen with ~2000 points per core, the expression gets compiled with numpy,
    significantly increasing computations times, while being compatible with the
    full computation precision up to 10e-9.
    """
    n_cores = mp.cpu_count()
    N = np.array(var_data).shape[-1]
    if N<(2000*n_cores):
        # Multicore computation
        print(f"Computation using all {n_cores} cores.")
        val_dicts = _data_dict_creation(variables, var_data)
        computed_func = _multicore_computation(n_cores, func, val_dicts)
    else:
        # lambdify computation
        print(
            "WARNING: computation time exceeding 30s. Compiling expression with NumPy."
            )
        computed_func = _lambdified_computation(func, variables, var_data)
    return np.array(computed_func)

def compute_error(err_func, variables, var_data, var_errors, corr_values:list=None):
    """
    Numerical computation of the error-formula for the input function.

    Parameters
    ----------
    func : sympy function
        The symbolic function for which the error needs to be computed.
    variables : list of sympy symbols
        The list of symbolic variables in the function.
    var_data : list of ndarray
        The list of data arrays corresponding to each variable.
    var_errors : list of ndarray
        The list of error arrays corresponding to each variable.
    corr_values : list, optional
        The list of correlation values between the variables. The default is None.

    Returns
    -------
    computed_error : list of floats
        The computed numerical error for the input function.

    """
    data = var_data+var_errors
    if corr_values is not None:
        data += corr_values
    computed_error = compute_numerical_function(err_func, variables, data)
    return computed_error

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
    try:
        assert len(corrs) % 2 == 0
        assert len(corrs) != 0
        for i in range(0,len(corrs), 2):
            eformula = sp.nsimplify(error_formula.subs(corrs[i+1], corrs[i]))
            corrs.pop(i+1)
    except AssertionError:
        eformula = error_formula
    erp = {
        "error_formula": eformula,
        "error_variables": {
            'variables': variables,
            'errors': errors
            }
    }
    if correlation:
        erp['error_variables']["corrs"] = corrs
    return erp

def gaus_legendre_integrator(fnc, a, b, points):
    """
    Integrates the function fcn(x) between a and b using the Gauss-Legendre method.

    Parameters
    ----------
    fcn : callable
        The function to integrate, written as a callable python function.
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
            - profiles: (normalized) density and w0 profiles with respect to the
            dimentionless radial distance from the centre
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
    result = subprocess.run([_king_exe, w0], capture_output=True, text=True,
                            cwd=_king_dir, check=False)
    if result.returncode != 0:
        print("Error during trhe Fortran90 code execution:")
        print(result.stderr)
    else:
        print("Calling Fortran90 code executor...\n")
        print(result.stdout)
    filelist = get_file_list(fold=_king_dir, key='.dat')
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

def _multicore_computation(n_cores, func, val_dicts):
    """
    Computation of the input function using multicore parallelization

    Parameters
    ----------
    n_cores : int
        Number of cores used.
    func : sympy function
        The function to numerically compute.
    val_dicts : dict
        The dictionaries of variable-values used to compute the function.

    Returns
    -------
    computed_func : ndarray
        Array of the computed function.

    """
    compute = __compute_sympy(func)
    with mp.Pool(n_cores) as pool:
        computed_func = pool.map(compute.compute, val_dicts)
    return computed_func

def _lambdified_computation(func, variables, var_data):
    """
    Compute the input function compiling the sympy expression using numpy, for
    way faster computation.

    Parameters
    ----------
    func : sympy function
        The function to compute.
    variables : list of sympy symbols
        The variables of the function.
    var_data : ndarray
        The values for the variables of the function.

    Returns
    -------
    computed_func : ndarray
        The result computed function.
    """
    f_lambdified = sp.lambdify(variables, func, modules="numpy")
    if not isinstance(var_data, list):
        var_data = [var_data]
    result = f_lambdified(*var_data)
    computed_func = []
    for x in result:
        computed_func.append(float(x))
    return computed_func

def _data_dict_creation(variables, var_data:ArrayLike):
    """
    function which creates the list of dictionaries in the format needed to compute
    sympy functions.

    Parameters
    ----------
    variables : list of sympy symbols
        The variables of the function to compute.
    var_data : ndarray
        the values for the variables of the function to compute.

    Returns
    -------
    val_dicts : list of dict
        A list of dictionaries, each containing the variables with an associated
        value, to be passed to the function in order to be computed.
    """
    N = np.array(var_data).shape[-1]
    val_dicts = []
    for n in range(N):
        data_dict = {}
        for i, var in enumerate(variables):
            var_name = f"{var}"
            data_dict[var_name] = var_data[i][n]
        val_dicts.append(data_dict)
    return val_dicts

class __compute_sympy():
    """
    Sub-Class for multiprocessing computation
    """
    def __init__(self, func):
        self.f = func
    def compute(self, vals):
        return float(sp.N(self.f.subs(vals)))
