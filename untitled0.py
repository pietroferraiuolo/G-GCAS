from ggcas import plots as gplt
import numpy as np
x = np.random.randn(1000)
z = np.random.randn(1000)
gplt.scat_xhist(x, z)
gplt.scat_xhist(x, z, xerr=0.051, xlabel=r'$v_{los}$', ylabel=r'$R$', color='red')