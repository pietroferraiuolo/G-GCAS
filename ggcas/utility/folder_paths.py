"""
Author(s) 
---------
    - Pietro Ferraiuolo : written in 2024

Description
-----------
"""
import os

try:
    BASE_DATA_PATH = os.environ['GCASDATA']
except KeyError as exc:
    raise KeyError("Environment variable not found! Define the GCASDATA\
                   env variable that points to '.../G-GCAS/ggcas/data") from exc
CATALOG_FILE            = os.path.join(BASE_DATA_PATH, 'Catalogue.xlsx')
QUERY_DATA_FOLDER       = os.path.join(BASE_DATA_PATH, 'query')
KING_MODELS_FOLDER      = os.path.join(BASE_DATA_PATH, 'king_model')

def CLUSTER_DATA_FOLDER(name:str):
    """
    Returns the cluster's data path.
    """
    path = os.path.join(QUERY_DATA_FOLDER, name.upper())
    return path

def CLUSTER_MODEL_FOLDER(name:str):
    """
    Returns the cluster's data path.
    """
    path = os.path.join(KING_MODELS_FOLDER, name.upper())
    return path