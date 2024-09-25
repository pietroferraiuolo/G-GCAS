from typing import List, Union, Dict, Any, Callable
import numpy as np
import sympy as sp
from astropy.table import QTable

def compute_numerical_function(
    func: sp.Basic, variables: List[sp.Symbol], var_data: List[np.ndarray]
) -> np.ndarray:
    """
    Compute the numerical value of a function, passing by the function, its
    variables and the data associated to the variables.
    """
    ...

def compute_error(
    func: sp.Basic,
    variables: List[sp.Symbol],
    var_data: List[np.ndarray],
    var_errors: List[np.ndarray],
    corr: bool = False,
    corr_values: List[np.ndarray] = None,
) -> np.ndarray:
    """
    Numerical computation of the error-formula for the input function.
    """
    ...

def error_propagation(
    func: sp.Basic, variables: List[sp.Symbol], correlation: bool = False
) -> Dict[str, Any]:
    """
    Computes the input function's error with the standard error propagation
    formula.
    """
    ...

def gaus_legendre_integrator(
    fnc: Callable[[float], float], a: float, b: float, points: int
) -> float:
    """
    Integrates the function fnc(x) between a and b using the Gauss-Legendre method.
    """
    ...

def king_integrator(
    w0: Union[float, int], output: str = "profile"
) -> Union[str, List[str]]:
    """
    This function calls a Fortran90 code for the Single-Mass King model integration
    routine.
    """
    ...

def _multicore_computation(
    n_cores: int, func: sp.Basic, val_dicts: List[Dict[str, float]]
) -> np.ndarray:
    """
    Computation of the input function using multicore parallelization.
    """
    ...

def _lambdified_computation(
    func: sp.Basic, variables: List[sp.Symbol], var_data: List[np.ndarray]
) -> np.ndarray:
    """
    Compute the input function compiling the sympy expression using numpy, for
    way faster computation.
    """
    ...

def _data_dict_creation(
    variables: List[sp.Symbol], var_data: List[np.ndarray]
) -> List[Dict[str, float]]:
    """
    Function which creates the list of dictionaries in the format needed to compute
    sympy functions.
    """
    ...

class __compute_sympy:
    """
    Support-Class for multiprocessing computation.
    """

    def __init__(self, func: sp.Basic) -> None: ...
    def compute(self, vals: Dict[str, float]) -> float: ...
