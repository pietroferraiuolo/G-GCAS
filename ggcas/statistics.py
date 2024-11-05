"""
Author(s)
---------
    - Pietro Ferraiuolo : Written in 2024

Description
-----------

"""
import os
import time
import numpy as np
import pandas as pd
from ggcas.analyzers._Rcode import check_packages
from ggcas._utility import R_SOURCE_FOLDER as _RSF
import rpy2.robjects as ro
from rpy2.robjects import (
    pandas2ri as pd2r,
    numpy2ri as np2r,
    r as R,
    globalenv as genv
)

def gaussian_mixture_estimator(data, n_components=1, verbose=False):
    """
    Gaussian Mixture Estimation function.
    """
    check_packages("mclust")
    np2r.activate()
    regression_code=os.path.join(_RSF, 'regression.R')
    R(f'source("{regression_code}")')
    reg_func = genv["regression"]
    r_data = np2r.numpy2rpy(data)
    r_result = reg_func(r_data, method='mclust', verb=verbose)
    x_gmm = np.array(r_result.rx2('x'))
    y_gmm = np.array(r_result.rx2('y'))
    coeffs = np.array(r_result.rx2('coeffs'))
    return x_gmm, y_gmm, coeffs

def kde_estimator(data, kind='gaussian', verbose=False):
    """
    Kernel Density Estimation function.
    """
    check_packages("minpack.lm")
    np2r.activate()
    regression_code=os.path.join(_RSF, 'regression.R')
    R(f'source("{regression_code}")')
    reg_func = genv["regression"]
    r_data = np2r.numpy2rpy(data)
    r_result = reg_func(r_data, method=kind, verb=verbose)
    x_kde = np.array(r_result.rx2('x'))
    y_kde = np.array(r_result.rx2('y'))
    coeffs = np.array(r_result.rx2('coeffs'))
    return x_kde, y_kde, coeffs

def _seed():
    return int(time.time())