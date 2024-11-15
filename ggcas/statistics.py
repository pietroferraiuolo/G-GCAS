"""
Author(s)
---------
    - Pietro Ferraiuolo : Written in 2024

Description
-----------

"""

import os as _os
import time as _time
import numpy as _np
import pandas as _pd
from ggcas.analyzers._Rcode import check_packages
from astropy.table import Table as _Table
from ggcas._utility import R_SOURCE_FOLDER as _RSF
import rpy2.robjects as _ro
from astroML.density_estimation import XDGMM
from rpy2.robjects import pandas2ri as pd2r, numpy2ri as np2r, r as R, globalenv as genv


def _seed():
    return int(_time.time())


def XD_estimator(data, errors, correlations=None, *xdargs):
    """
    Extreme Deconvolution Estimation function.

    Parameters
    ----------
    data : numpy.ndarray
        The data to be analyzed.
    errors : numpy.ndarray
        The errors associated with the data.
    correlations : dict
        The correlations between the data. If a dictionary, the keys are
        (i, j) tuples representing feature pairs and the values are
        (n_samples,) arrays of correlation coefficients for each sample.
    *xdargs : optional
        XDGMM hyper-parameters for model tuning. See
        <a href="https://www.astroml.org/modules/generated/astroML.density_estimation.XDGMM.html">
        astroML</a> documentation for more information.

    Returns
    -------
    model : astroML.density_estimation.XDGMM
        The XDGMM fitted model.
    """
    x_data = _data_format_check(data)  # (n_samples, n_features)
    errors = _data_format_check(errors)  # (n_samples, n_features, n_features)
    if correlations is not None:
        covariance_matrix = _construct_covariance_matrices(errors, correlations)
    else:
        covariance_matrix = _np.array([_np.diag(e**2) for e in errors])
    model = XDGMM(random_state=_seed(), *xdargs)
    model.fit(x_data, covariance_matrix)
    return model


def gaussian_mixture_model(train_data, fit_data=None, **kwargs):
    """
    Gaussian Mixture Estimation function.

    Parameters
    ----------
    data : numpy.ndarray
        The data to be fitted with a gaussian mixture model.
    **kwargs : optional
        Additional keyword arguments for the R GM_model function. See
        <a href="https://cran.r-project.org/web/packages/mclust/mclust.pdf">
        mclust</a> documentation for more information.

    Returns
    -------
    fitted_model : dict
        The fitted gaussian mixture model and its parameters.
    """
    check_packages("mclust")
    np2r.activate()
    code = _os.path.join(_RSF, "gaussian_mixture.R")
    R(f'source("{code}")')
    if fit_data is not None:
        r_data = np2r.numpy2rpy(train_data)
        r_fit_data = np2r.numpy2rpy(fit_data)
        # Convert kwargs to R list
        r_kwargs = _ro.vectors.ListVector(kwargs)
        # Call the R function with the data and additional arguments
        fitted_model = genv["GaussianMixtureModel"](r_data, r_fit_data, **dict(r_kwargs.items()))
        clusters = fitted_model.rx2("cluster")
        fitted_model = fitted_model.rx2("model")
    else:
        r_data = np2r.numpy2rpy(train_data)
        # Convert kwargs to R list
        r_kwargs = _ro.vectors.ListVector(kwargs)
        # Call the R function with the data and additional arguments
        fitted_model = genv["GM_model"](r_data, **dict(r_kwargs.items()))
        clusters = None
    return _RGMModel(fitted_model, clusters)


def kde_estimator(data, kind="gaussian", verbose=False):
    """
    Kernel Density Estimation function.
    """
    check_packages("minpack.lm")
    np2r.activate()
    regression_code = _os.path.join(_RSF, "regression.R")
    R(f'source("{regression_code}")')
    reg_func = genv["regression"]
    r_data = np2r.numpy2rpy(data)
    r_result = reg_func(r_data, method=kind, verb=verbose)
    x_kde = _np.array(r_result.rx2("x"))
    y_kde = _np.array(r_result.rx2("y"))
    coeffs = _np.array(r_result.rx2("coeffs"))
    return x_kde, y_kde, coeffs


def _data_format_check(data):
    """
    Function which checks and formats the input data to be ready
    for the XDGMM model fit.

    Parameters:
    ----------
    data : numpy.ndarray, list, astropy.table.Table, pandas.DataFrame
        The data whose format has to be checked.
    
    Returns:
    -------
    data : numpy.ndarray
        The data in the correct format for the XDGMM model.
        Rturned in shape (n_samples, n_features).
    """
    if isinstance(data, (_Table, _pd.DataFrame)):
        data = data.to_numpy()
    elif isinstance(data, list):
        data = _np.stack(data).T
    return data


def _construct_covariance_matrices(errors, correlations: dict):
    """
    Constructs covariance matrices for each sample based on given errors and correlations.

    Parameters:
    ----------
    errors : numpy.ndarray
        An (n_samples, n_features) array with standard deviations for each feature.
    correlations : dict
        A dictionary of correlation arrays where keys are (i, j) tuples
        representing feature pairs and values are (n_samples,) arrays of
        correlation coefficients for each sample.

    Returns:
    -------
    cov_mat : numpy.ndarray
        Array of covariance matrices, one for each sample.
    """
    if not isinstance(correlations, dict):
        raise ValueError(
            f"Correlations must be a dictionary, not {type(correlations)}."
        )
    n_samples, n_features = errors.shape
    X_error = []
    for i in range(n_samples):
        # Initialize covariance matrix for sample i
        cov_matrix = _np.zeros((n_features, n_features))
        # Set variances on the diagonal
        for j in range(n_features):
            cov_matrix[j, j] = errors[i, j] ** 2
        # Set covariances for each pair based on correlations
        for (f1, f2), rho in correlations.items():
            cov_matrix[f1, f2] = rho[i] * errors[i, f1] * errors[i, f2]
            cov_matrix[f2, f1] = cov_matrix[f1, f2]  # symmetry
        X_error.append(cov_matrix)
    cov_mat = _np.array(X_error)
    return cov_mat


class _RGMModel:
    """
    Class to convert an R model to a Python dictionary.
    """

    def __init__(self, r_model, predictions=None):
        """The Constructor"""
        self.rmodel         = r_model
        self.model          = self._listvector_to_dict(r_model)
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
        return _np.array([self.model["n"], self.model["d"], self.model["G"]])

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
        return _np.array([self.model["z"], self.model["classification"]])
    

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
        self.model["parameters"]["variance"]["covMats"] = _np.swapaxes(
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

    @staticmethod
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
                    _ro.rinterface.IntSexpVector,
                    _ro.rinterface.FloatSexpVector,
                    _ro.rinterface.StrSexpVector,
                ),
            ):
                py_dict[key] = list(value) if len(value) > 1 else value[0]
            # Handle data frames (convert to pandas DataFrame)
            elif isinstance(value, _ro.vectors.DataFrame):
                py_dict[key] = pd2r.rpy2py(value)
            # Handle nested ListVectors
            elif isinstance(value, _ro.vectors.ListVector):
                py_dict[key] = _RGMModel._listvector_to_dict(value)
            # Handle other lists
            elif isinstance(value, _ro.vectors.Vector):
                py_dict[key] = list(value)
            # Other R objects could be added as necessary
            else:
                py_dict[key] = value
        return py_dict
