"""
Author(s)
---------
    - Pietro Ferraiuolo : Written in 2024

Description
-----------
A module which contains the Cluster class, which contains all the information of
a specified cluster.

How to Use
----------
Initialize the class with a cluster's name. As example

>>> from ggcas.cluster import Cluster
>>> ngc104 = Cluster('ngc104')

Now we can call methods to read the parameters

>>> ngc104.id
'NGC104'
>>> ngc104.w0
8.82
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.table import Table
from ggcas.utility import folder_paths as fn

class Cluster:
    """
    Class for the cluster parameter loading.

    Upon initializing, it is loaded with the specified cluster's parameters, ta
    ken from the Harris Catalogue 2010 Edition.

    Methods
    -------
    _loadClusterParameters(self, name) : function
        Loads the desired cluster's parameters

    How to Use
    ----------
    Initialize the class with a cluster's name. As example

    >>> from ggcas.cluster import Cluster
    >>> ngc104 = Cluster('ngc104')
    """
    def __init__(self, name: str):
        """The constructor"""
        self.id     = name.upper()
        self.data_path = fn.CLUSTER_DATA_FOLDER(self.id)
        parms       = self._load_cluster_parameters(self.id)
        self.model  = self._load_king_model()
        self.ra     = parms.loc['ra']*u.deg
        self.dec    = parms.loc['dec']*u.deg
        self.dist   = parms.loc['dist']*u.kpc
        self.rc     = parms.loc['rc']/60 * u.deg
        self.rh     = parms.loc['rh']/60 * u.deg
        self.w0     = parms.loc['w0']
        self.logc   = parms.loc['logc']
        self.cflag  = parms.loc['collapsed']=='Y'

    def show_model(self, **kwargs):
        """
        Function for plotting the loaded king model.

        Optional Parameters
        -------------------
        **kwargs :
            color : color of the main plot.
            scale : scale of the axes, default linear, can be set to 'log'.
        """
        scale = kwargs.get('scale', None)
        c = kwargs.get('color', 'black')
        plt.figure(figsize=(8,6))
        plt.plot(self.model['xi'], self.model['w'], color=c)
        plt.plot([self.model['xi'].min(), self.model['xi'].min()],
                 [self.model['w'].min(), self.model['w'].max()],
                 c='red', linestyle='--',
                 label=rf"$W_0$={self.model['w'].max()}")
        plt.xlabel(r"$\xi$ = $\dfrac{r}{r_t}$")
        plt.ylabel("w")
        plt.title('Integrated King Model')
        if scale=='log':
            plt.xscale('log')
            plt.yscale('log')
        plt.legend(loc='best')
        plt.show()

    def _load_cluster_parameters(self, name: str):
        """
        Loads the parameters of the requested cluster from the Harris Catalog
        2010 Globular Cluster Database, written in the Catalogue.xlsx file

        Parameters
        ----------
        name : str
            Name of the requested Globular Cluster.

        Returns
        -------
        cat_row : TYPE
            Pandas Series with all the necessary paramenters to ilitialize the Cluster Class.
        """
        catalog = pd.read_excel(fn.CATALOG_FILE, index_col=0)
        cat_row = catalog.loc[name.upper()]
        return cat_row

    def _load_king_model(self):
        """
        Loads the integrated model of the cluster.

        Parameters
        ----------
        kind : str, optional
            Model to load. The default is 'SM', that stands for 'Single Mass'.
            Values are:
                SM: Single Mass King model.

        Returns
        -------
        model : astropy table
            Astropy table containing the integrated quantities of the king model.
            These are:
                'xi': dimentionless radial distance from gc center, normalized
                at tidal radius.
                'w': the w-king parameter, that is essentially the gravitational
                potential.
                'rho': the noralized density profile of the clusted.
        """
        model = Table()
        file = os.path.join(fn.CLUSTER_MODEL_FOLDER(self.id),'SM_king.txt')
        model['xi']     = np.loadtxt(file, skiprows=1, usecols=1)
        model['w']      = np.loadtxt(file, skiprows=1, usecols=2)
        model['rho']    = np.loadtxt(file, skiprows=1, usecols=3)
        return model
