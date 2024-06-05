"""
Created on May 2024
    -Author: P.Ferraiuolo
"""
import sympy as sp
import numpy as np

def errPropagation(func, variables, vars_values,  compute=True):
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

    if compute==True:
        values = dict()
        for xx in range(len(variables)):
            values[variables[xx]] = 0

        computed_error = np.zeros(len(vars_values[0]))
        for ii in range(len(vars_values[0])):
            for xx in range(len(variables)):
                values[variables[xx]] = vars_values[xx][ii]

            computed_error[ii] = sp.N(error_formula.subs(values))

        return computed_error
    else:
        return error_formula

