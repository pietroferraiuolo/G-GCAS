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
import rpy2.robjects as ro
from rpy2.robjects import (
    pandas2ri as pd2r,
    numpy2ri as np2r,
    r as R,
    globalenv as genv
)
from ggcas._utility import R_SOURCE_FODLER as _RSF

def _seed():
    return int(time.time())