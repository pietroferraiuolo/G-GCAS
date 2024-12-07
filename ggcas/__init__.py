"""
GAIA - GLOBULAR CLUSTER ANALYSIS SOFTWARE
=========================================
Author(s)
---------
    - Pietro Ferraiuolo : Written in 2024

Description
-----------

"""
from . import statistics as stat
from . import plots as plt
from ._query import (
    GaiaQuery, 
    available_tables
)
from ._cluster import Cluster
from ._utility import (
    load_data,
    get_file_list,
    tnlist
)
from .analyzers import calculus as calc
from .analyzers import dynamics as dyn
from .analyzers import mcluster as mcl

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
