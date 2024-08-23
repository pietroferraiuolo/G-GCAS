"""
Author(s) 
---------
    - Pietro Ferraiuolo : written in 2024

Description
-----------
"""
import time
runfile('/home/pietrof/git/G-GCAS/ggcas/initGGCAS.py', wdir='/home/pietrof/git/G-GCAS/ggcas')
n6121 = Cluster('ngc6121')
cond = ['parallax IS NOT NULL', 'pm IS NOT NULL']
# m4_p = dr3.get_photometry(radius=n6121.rt, gc=n6121, conds=cond)
m4_a = dr3.get_astrometry(radius=n6121.rt, gc=n6121, conds=cond)
# hist = gplt.histogram(m4_a.sample['parallax'], kde=True, xlim=(-5,5))
w = gfnc.angular_separation(n6121.ra, n6121.dec)
var_data = [m4_a.sample['ra'], m4_a.sample['dec']]
start = time.time()
computed_w = calculus.compute_numerical_function(w['f'], w['vars'], var_data)
stop = time.time()
ct = stop-start
print(f'Computation time: {ct:.3f} s -> {ct/60:.2f} min')
