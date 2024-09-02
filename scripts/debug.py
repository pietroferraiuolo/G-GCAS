"""
Author(s)
---------
    - Pietro Ferraiuolo : written in 2024

Description
-----------
"""
from time import time as t
import os, numpy as np
from matplotlib import pyplot as plt
from ggcas.cluster import Cluster
from ggcas.analyzers import calculus as gcalc
from ggcas.analyzers import effective_potential as gep
from ggcas.analyzers import dynamics as gdy
from ggcas.query import GaiaQuery
from ggcas import plots as gplt, functions as gfnc
dr3 = GaiaQuery()

def et(start, stop):
    print(f'{stop-start:.3f} s')

# =============================================================================
# Angular separation computation
# =============================================================================
def compute_angular_separation(sample):
    w = gfnc.angular_separation(n6121.ra, n6121.dec)
    data = [sample.sample['ra'], sample.sample['dec']]
    errors = [sample.sample['ra_error'], sample.sample['dec_error']]
    s = t()
    n_w = gcalc.compute_numerical_function(w['f'], w['vars'], data)
    ss = t()
    et(s,ss)
    s = t()
    c_ew = gcalc.compute_error(w['f'], w['vars'], data, errors)
    ss = t()
    et(s,ss)
    values = []
    errors = []
    for i, w in enumerate(n_w):
        if c_ew[i]/w<2:
            values.append(w)
            errors.append(c_ew[i])
    W = np.array([values,errors])
    return W
#______________________________________________________________________________
n6121 = Cluster('ngc6121')
cond = ['parallax IS NOT NULL', 'pm IS NOT NULL']
params = ['source_id', 'parallax', 'parallax_error', 'ra', 'ra_error', 'dec', 'dec_error']
m4_a = dr3.free_query(radius=n6121.rt, gc=n6121, data=params, conds=cond)

w = compute_angular_separation(m4_a)
wh = gplt.histogram(w[1]/w[0]*100, kde=True, xlabel='% error')
