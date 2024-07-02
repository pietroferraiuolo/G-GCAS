import os
from ggcas.cluster import Cluster
from ggcas import utils as gutils
from ggcas import plots as gplt
from ggcas import functions as gfunc
from ggcas.analyzer import *
from ggcas.query import GaiaQuery 

gq = GaiaQuery()

print("""

GAIA - GLOBULAR CLUSTERS ANALYSIS SOFTWARE
            ..............                                                                         
         ..:;;..:;;;;;:::::;;                                                                      
       ;;;;;;::.::;;;;;;;;;;;;;                                                                    
      ;;;;.:;;;..;XXXXXX.::....:                                                                   
     :;::;::..+XXXXXXXXX+:;;;;;;:                                                   
    ;::;:.:;;:XXXXXXXXXXX::::::::;        __ _  __ _(_) __ _                                       
    .;;..;;;:.:XXXXXXXXX$$$$$$$$$$X.     / _` |/ _` | |/ _` |                                      
    :;..:;;;..:xXXXXXXX$$$$$$$$$$$$X    | (_| | (_| | | (_| |                                      
    :;:.:;;.XXXXXXXXX$$$$$$$$$$$$$$;     \__, |\__,_|_|\__,_|                                      
    .;;:.:X$$$$$$$$$$$$$$$$$$$$$$X.      |___/                                      
    ..:;:$$$$$$$$$$$$$$$$$$$$$$X;.                                                                 
     :;;;$$$$$$$$$$$$$$$$$$$$::;;                                                                  
      ...;$$$$$$$$$$$$$$x;:;;;;;                                                                   
        ......:;:....;;;;;;;:.                                                                     
         ::::::::::;;;::...                                                                       

             INITIALIZED                                                        

""")

try:
    print('Data path is {}'.format(os.environ['PYGCASCONF']))
except KeyError:
    raise KeyError("Environment variable not found. Please set the PYGCASCONF env variable that points to the data folder: '.../G-Gcas/ggcas/data/'")

def help():
    TEXT="""
G-GCAS PACKAGE DOCUMENTATION

- cluster

- functions

- query

- utils

- plots

- analysis
    """
    
    print(TEXT)
