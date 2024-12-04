"""
Author(s)
---------
- Pietro Ferraiuolo : Written in 2024

Description
-----------
This module containg a series of classes which are used to convert R models
into Python dictionaries. This allows for a easier and more intuitive access
to the model parameters and results.

"""
import numpy as np
import rpy2.robjects as ro
from rpy2.robjects import (
    pandas2ri as pd2r, 
    numpy2ri as np2r, 
    r as R,
)

class GMModel:
    """
    Class to convert the R Mclust Gaussian Mixture Model into a Python 
    dictionary.
    """

    def __init__(self, r_model, predictions=None):
        """The Constructor"""
        self.rmodel         = r_model
        self.model          = _listvector_to_dict(r_model)
        self.classification = predictions
        self._predicted     = False if predictions is None else True

    def __str__(self):
        str_ = \
f"""----------------------------------------------------
Gaussian finite mixture model fitted by EM algorithm
----------------------------------------------------

Mclust {self.model['modelName']} model fitted with {self.ndG[2]} components:

log-likelihood : {self.loglik}
n : {self.ndG[0]}
df : {self.ndG[1]}
BIC : {self.bic['BIC']}
ICL : {self.model['icl']}

Predicted : {self._predicted}
"""
        return str_
    
    def __repr__(self):
        repr_ = \
f"""Python wrapper for R Mclust Gaussian Mixture Model
--------------------------------------------------
.rmodel : R model object as rpy2.robjects
.model  : R Model translation into py dict

"""
        return repr_

    @property
    def data(self):
        """
        The data used to fit the model.
        """
        return self.model["data"]

    @property
    def bic(self):
        """
        The Bayesian Information Criterion of the model.
        """
        return {"bic": self.model["bic"], "BIC": self.model["BIC"]}

    @property
    def ndG(self):
        """
        Array containing in order:
            n: the number of data points used to fit the model
            d: the number of parameters/features fitted with the model
            G: the number of components in the model
        """
        return np.array([self.model["n"], self.model["d"], self.model["G"]])

    @property
    def loglik(self):
        """
        The log-likelihood of the model.
        """
        return self.model["loglik"]

    @property
    def train_classification(self):
        """
        Array containing in order:
            z: the membership probability of each data point to each component
            classification: the classification of the data points
        """
        return np.array([self.model["z"], self.model["classification"]])

    @property
    def parameters(self):
        """
        Dictionary containint all the parameters estimated by the model, for each component
        Components:
            pro: the mixing proportions
            mean: the means
            variance: another dictionary containing the variance and covariance matrices:
                covMats: the covariance matrices for each group of the model
                sigmasq: the variance of the features of the model
        """
        self.model["parameters"]["variance"]["covMats"] = np.swapaxes(
            self.model["parameters"]["variance"]["sigma"], 0, 2
        )
        keys_to_remove = ["modelName", "d", "G", "sigma", "scale"]
        for key in keys_to_remove:
            self.model["parameters"].pop(key, None)
        return self.model["parameters"]

    @property
    def uncertainty(self):
        """
        The uncertainty of the model for each data point membership probability.
        """
        return self.model["uncertainty"]


class RegressionModel:
    """
    Class to convert the R LM Regression Model into a Python dictionary.
    """

    def __init__(self, r_model, type:str='Unknown Model'):
        """The Constructor"""
        self.rmodel = r_model
        self.model  = _listvector_to_dict(r_model)
        self._model_kind = type


    def __repr__(self):
        """The representation of the model."""
        repr_ = \
f"""----------------------------------------------------
Python wrapper for R Levenberg-Marquardt Nonlinear 
              Least-Squares Algorithm            
----------------------------------------------------
.rmodel : R model object as rpy2.robjects
.model  : R Model translation into py dict

"""
        return repr_
    
    def __str__(self):
        """The string representation of the model."""
        txt = _kde_labels(self._model_kind, self.coeffs).splitlines()
        txt.pop(0)
        txt = '\n'.join(txt)
        str_ = \
f"""----------------------------------------------------------
{self._model_kind.upper()} Regression Model fitted by LM algorithm
----------------------------------------------------------
{txt}
"""
        return str_
    
    @property
    def coeffs(self):
        """
        The coefficients of the regression model.
        """
        return self.model["coeffs"]
    
    @property
    def data(self):
        """
        The data used to fit the model.
        """
        return self.model["data"]
    
    @property
    def x(self):
        """
        The independent variables of the model.
        """
        return self.model["x"]

    @property
    def y(self):
        """
        The dependent variables of the model.
        """
        return self.model["y"]
    
    @property
    def residuals(self):
        """
        The residuals of the model.
        """
        return self.model["residuals"]

    @property
    def kind(self):
        """
        The kind of model fitted to the data.
        """
        return self._model_kind


def _listvector_to_dict(r_listvector):
    """
    Recursively converts an R ListVector (from rpy2) to a nested Python dictionary.
    """
    py_dict = {}
    for key, value in r_listvector.items():
        # Handle simple types
        if isinstance(
            value,
            (
                ro.rinterface.IntSexpVector,
                ro.rinterface.FloatSexpVector,
                ro.rinterface.StrSexpVector,
            ),
        ):
            py_dict[key] = list(value) if len(value) > 1 else value[0]
        # Handle data frames (convert to pandas DataFrame)
        elif isinstance(value, ro.vectors.DataFrame):
            py_dict[key] = pd2r.rpy2py(value)
        # Handle nested ListVectors
        elif isinstance(value, ro.vectors.ListVector):
            py_dict[key] = _listvector_to_dict(value)
        # Handle other lists
        elif isinstance(value, ro.vectors.Vector):
            py_dict[key] = list(value)
        # Other R objects could be added as necessary
        else:
            py_dict[key] = value
    return py_dict

def _kde_labels(kind: str, coeffs):
    """
    Return the labels for the KDE plot.
    """
    def _format_number(num):
        """
        Format the number using scientific notation if it is too large or too small.
        """
        if abs(num) < 1e-3 or abs(num) > 1e3:
            return f"{num:.2e}"
        else:
            return f"{num:.2f}"
    
    if kind == 'gaussian':
        A, mu, sigma2 = coeffs
        label = f"""Gaussian
$A$   = {_format_number(A)}
$\mu$   = {_format_number(mu)}
$\sigma^2$  = {_format_number(sigma2)}"""
        
    elif kind == 'boltzmann':
        A1, A2, x0, dx = coeffs
        label = f"""Boltzmann
$A1$   = {_format_number(A1)}
$A2$   = {_format_number(A2)}
$x_0$   = {_format_number(x0)}
$dx$   = {_format_number(dx)}"""
        
    elif kind == 'exponential':
        A, lmbda = coeffs
        label = f"""Exponential
$A$   = {_format_number(A)}
$\lambda$ = {_format_number(lmbda)}"""
        
    elif kind == 'king':
        A, ve, sigma = coeffs
        label = f"""King
$A$   = {_format_number(A)}
$v_e$   = {_format_number(ve)}
$\sigma$  = {_format_number(sigma)}"""
        
    elif kind == 'maxwell':
        A, sigma = coeffs
        label = f"""Maxwell
$A$   = {_format_number(A)}
$\sigma$  = {_format_number(sigma)}"""
        
    elif kind == 'rayleigh':
        A, sigma = coeffs
        label = f"""Rayleigh
$A$   = {_format_number(A)}
$\sigma$  = {_format_number(sigma)}"""
        
    elif kind == 'lorentzian':
        A, x0, gamma = coeffs
        label = f"""Lorentzian
$A$   = {_format_number(A)}
$x_0$   = {_format_number(x0)}
$\gamma$  = {_format_number(gamma)}"""
    
    elif kind == 'power':
        A, n = coeffs
        label = f"""Power
$A$   = {_format_number(A)}
$n$   = {_format_number(n)}"""
        
    elif kind == 'lognormal':
        A, mu, sigma = coeffs
        label = f"""Lognormal
$A$   = {_format_number(A)}
$\mu$   = {_format_number(mu)}
$\sigma$  = {_format_number(sigma)}"""
    
    elif kind == 'poisson':
        A, lmbda = coeffs
        label = f"""Poisson
$A$   = {_format_number(A)}
$\lambda$ = {_format_number(lmbda)}"""
        
    elif kind == 'linear':
        A, B = coeffs
        label = f"""Linear
$A$   = {_format_number(A)}
$B$   = {_format_number(B)}"""
    
    return label
