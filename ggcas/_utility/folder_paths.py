"""
Author(s)
---------
    - Pietro Ferraiuolo : written in 2024

Description
-----------
"""
import os

try:
    BASE_PATH = os.getenv('GCASDATA')
except KeyError as exc:
    raise KeyError("Environment variable not found! Define the GCASDATA env variable that points to '.../G-GCAS/ggcas/data") from exc
BASE_DATA_PATH          = os.path.join(BASE_PATH, 'data')
CATALOG_FILE            = os.path.join(BASE_PATH, '_Catalogue.xlsx')
QUERY_DATA_FOLDER       = os.path.join(BASE_DATA_PATH, 'query')
KING_MODELS_FOLDER      = os.path.join(BASE_DATA_PATH, 'models')
KING_INTEGRATOR_FOLDER  = os.path.join(BASE_PATH, 'analyzers', '_king')
SIMULATION_FOLDER       = os.path.join(BASE_DATA_PATH, 'simulations')
R_SOURCE_FOLDER         = os.path.join(BASE_PATH, 'analyzers', '_Rcode')
MCLUSTER_SOURCE_CODE    = os.path.join(BASE_PATH, 'analyzers', '_mcluster')
UNTRACKED_DATA_FOLDER   = os.path.join(BASE_DATA_PATH,'UntrackedData')

paths = [BASE_DATA_PATH, QUERY_DATA_FOLDER, KING_MODELS_FOLDER, SIMULATION_FOLDER, UNTRACKED_DATA_FOLDER]
for p in paths:
    if not os.path.exists(p):
        os.mkdir(p)

def CLUSTER_DATA_FOLDER(name:str):
    """
    Returns the cluster's data path.
    """
    path = os.path.join(QUERY_DATA_FOLDER, name.upper())
    return path

def CLUSTER_MODEL_FOLDER(name:str):
    """
    Returns the cluster's model path.
    """
    path = os.path.join(KING_MODELS_FOLDER, name.upper())
    return path
