"""
Author(s)
---------
    - Pietro Ferraiuolo : written in 2024

Description
-----------
This module contains R libraries implementation checks
for the G-GCAS package.
"""
from typing import Union
from rpy2.robjects.packages import importr
from rpy2.rinterface_lib.callbacks import logger as rpy2_logger
import logging
# Disable unnecessary R logging
rpy2_logger.setLevel(logging.ERROR)

utils = importr('utils')
utils.chooseCRANmirror(ind=1)

def check_packages(packages:Union[str, list[str]]) -> None:
    """
    Check if the R packages are installed.

    Parameters
    ----------
    packages : str or list[str]
        The name of the R package(s) to check.

    Returns
    -------
    None
    """
    if isinstance(packages, str):
        packages = [packages]
    for package in packages:
        try:
            importr(package)
        except:
            print(f"Package '{package}' is not installed.\nInstalling it now...")
            utils.install_packages(package)
            importr(package)
            print(f"Package '{package}' installed.")