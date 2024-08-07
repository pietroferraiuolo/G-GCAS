"""
Author(s):
----------
    - Pietro Ferraiuolo: written in 2024

Description
-----------

How to Use it
-------------

Documentation
=============

"""
import os, numpy as np
from matplotlib import pyplot as plt
from ggcas.cluster import Cluster
from ggcas.analyzers import calculus, effective_potential as efp, kinematics as kin
from ggcas.query import GaiaQuery
from ggcas import plots as gplt, functions as gfnc
dr3 = GaiaQuery()

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
