import os
from ggcas.cluster import Cluster
from ggcas import plots as gplt
#from ggcas import funcions, analyzer
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

print('Data path is {}'.format(os.environ['PYGCASCONF']))
<<<<<<< Updated upstream
=======

def help():
    TEXT="""
G-GCAS PACKAGE DOCUMENTATION

- cluster

- functions

- query

- plots

- analysis
    """
    
    print(TEXT)
    
>>>>>>> Stashed changes
