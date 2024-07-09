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
import astropy.units as u
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
        parms       = self._load_cluster_parameters(name.upper())
        self.king   = self._load_king_model()
        self.ra     = parms.loc['ra']*u.deg
        self.dec    = parms.loc['dec']*u.deg
        self.dist   = parms.loc['dist']*u.kpc
        self.rc     = parms.loc['rc']/60 * u.deg
        self.rh     = parms.loc['rh']/60 * u.deg
        self.w0     = parms.loc['w0']
        self.logc   = parms.loc['logc']
        self.cflag  = parms.loc['collapsed']=='Y'

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
        
        return