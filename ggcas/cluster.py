# -*- coding: utf-8 -*-
"""
Created on May 2024
    -Author: P.Ferraiuolo
"""
import os
import pandas as pd
import astropy.units as u

class Cluster:
    
    def __init__(self, name: str):
        '''The constructor'''
        parms       = self._loadClusterParameters(name.upper())
        self.id     = name.upper()
        self.ra     = parms.loc['ra']*u.deg
        self.dec    = parms.loc['dec']*u.deg
        self.dist   = parms.loc['dist']*u.kpc
        self.rc     = parms.loc['rc']/60 * u.deg
        self.rh     = parms.loc['rh']/60 * u.deg
        self.w0     = parms.loc['w0']
        self.logc   = parms.loc['logc']
        self.cflag  = parms.loc['collapsed']=='Y'
        
    
    def _loadClusterParameters(self, name: str):
        '''
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
        '''
        path = os.environ['PYGCASCONF']
        catalog = pd.read_excel(os.path.join(path, 'Catalogue.xlsx'), index_col=0)
        cat_row = catalog.loc[name.upper()]
        return cat_row
    
