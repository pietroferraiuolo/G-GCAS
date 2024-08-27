"""
Author(s) 
---------
    - Pietro Ferraiuolo : written in 2024

Description
-----------
"""
runfile('/home/pietrof/git/G-GCAS/ggcas/initGGCAS.py', wdir='/home/pietrof/git/G-GCAS/ggcas')
n6121 = Cluster('ngc6121')
cond = ['parallax IS NOT NULL', 'pm IS NOT NULL']
params = ['source_id', 'parallax', 'parallax_error', 'ra', 'ra_error', 'dec', 'dec_error']
# m4_p = dr3.get_photometry(radius=n6121.rt, gc=n6121, conds=cond)
m4_a = dr3.free_query(radius=n6121.rt, gc=n6121, data=params, conds=cond)
# hist = gplt.histogram(m4_a.sample['parallax'], kde=True, xlim=(-5,5))
w = gfnc.angular_separation(n6121.ra, n6121.dec)
data = [m4_a.sample['ra'], m4_a.sample['dec']]
errors = [m4_a.sample['ra_error'], m4_a.sample['dec_error']]
n_w = gcalc.compute_numerical_function(w['f'], w['vars'], data)
e_w = gcalc.error_propagation(w['f'], w['vars'])
c_ew = gcalc.compute_error(w['f'], w['vars'], data, errors)
W = np.array([n_w, c_ew])
values = []
errors = []
for i, w in enumerate(n_w):
    if c_ew[i]/w<2:
        values.append(w)
        errors.append(c_ew[i])
W = np.array([values,errors])
wh = gplt.histogram(W[1]/W[0]*100, kde=True, xlabel='% error')
