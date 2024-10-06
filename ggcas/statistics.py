"""
Author(s)
---------
    - Pietro Ferraiuolo : Written in 2024

Description
-----------

"""
import time
import numpy as np
import pandas as pd
import emcee
import pymc as pm
from astroML import XDGMM
from sklearn.mixture import (
    BayesianGaussianMixture,
    GaussianMixture
)
from sklearn.metrics import (
    classification_report as cr,
    confusion_matrix as cm,
    accuracy_score as acs
)
from sklearn.model_selection import train_test_split as tts

def gaussianMixture(data, *args, **kwargs):
    rs = _seed()
    gm = GaussianMixture(*args, **kwargs)
    train_data, temp_data = tts(data, train_size=0.9, random_state=rs)
    validation_data, test_data = tts(temp_data, test_size=0.5, random_state=rs)
    gm.fit(train_data)
    return gm, validation_data, test_data

def bayesianGaussianMixture(data, *args, **kwargs):
    rs = _seed()
    gm = BayesianGaussianMixture(*args, **kwargs)
    train_data, temp_data = tts(data, train_size=0.9, random_state=rs)
    validation_data, test_data = tts(temp_data, test_size=0.5, random_state=rs)
    gm.fit(train_data)
    return gm, validation_data, test_data

def _seed():
    return int(time.time())