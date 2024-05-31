# -*- coding: utf-8 -*-
"""
Created on May 2024
    -Author: P. Ferraiuolo
"""
import os
import datetime as dt
from astropy.table import Table
import astropy.units as u
from astroquery.gaia import Gaia
from typing import Optional, Union


class GaiaQuery:
    
    def __init__(self, gaia_table: Optional[Union[str, list]] = "gaiadr3.gaia_source"):
        '''The Constructor'''
        Gaia.MAIN_GAIA_TABLE = gaia_table
        Gaia.ROW_LIMIT = -1
        
        self._table = gaia_table
        self._q     = Gaia
        self._path  = "C:/Users/Er_da/Desktop/G-GCAS/data/query/"
        self._baseQ = """SELECT {data}
                FROM {table} 
                WHERE CONTAINS(POINT('ICRS',gaiadr3.gaia_source.ra,gaiadr3.gaia_source.dec),CIRCLE('ICRS',{circle}))=1
                  {cond}
                """
        self._joinQ = """SELECT {data}
                FROM {table} 
                WHERE CONTAINS(POINT('ICRS',gaiadr3.gaia_source.ra,gaiadr3.gaia_source.dec),CIRCLE('ICRS',{circle}))=1
                  {cond}
                """
        
    def _tn(self):
        '''
        Returns a tracking number, with format YYYYMMDD_HHMMSS

        Returns
        -------
        tn : str
            Tracking number as string.
        '''
        tn = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
        return tn
    
    def _checkPathExist(self, dest: str):
        '''
        

        Parameters
        ----------
        dest : TYPE
            DESCRIPTION.

        Returns
        -------
        path : TYPE
            DESCRIPTION.

        '''
        fold = self._path+dest.upper()+'/'
        if not os.path.exists(fold):
            os.makedirs(fold)
            print(f"Path '{fold}' did not exist. Created.")
            return fold
        else:
            return fold
        
    def _adqlWriter(self, ra, dec, radius, data: Optional[Union[str,list]] = None, conditions: Optional[Union[str,list]] = None):
        '''
        

        Parameters
        ----------
        ra : TYPE
            DESCRIPTION.
        dec : TYPE
            DESCRIPTION.
        radius : TYPE
            DESCRIPTION.
        data : Optional[Union[str,list]], optional
            DESCRIPTION. The default is None.
        conditions : Optional[Union[str,list]], optional
            DESCRIPTION. The default is None.

        Returns
        -------
        query : TYPE
            DESCRIPTION.

        '''
        
        if isinstance(ra, u.Quantity) and isinstance(dec, u.Quantity):
            ra = str(ra/u.deg)
            dec= str(dec/u.deg)
        else:
            ra = str(ra)
            dec = str(dec)
            
        radius  = str(radius)
        
        circle  = ra + "," + dec + "," + radius
        dat     = ''
        cond    = ''
        
        if data is not None:
            if isinstance(data, list):
                for i in range(len(data)-1):
                    dat += data[i]+', '
                dat += data[len(data)-1]
            else: dat=data
        
        if conditions is not None:
            if isinstance(conditions, str):
                conditions = conditions.split(',')
            cond = '    AND '
            for i in range(len(conditions)-1):
                cond += conditions[i]+"""
            AND """
            cond += conditions[len(conditions)-1]
        
        query = self._baseQ.format(data=dat, table=self._table, circle=circle, cond=cond)
        
        return query
        
    
    def freeQuery(self, ra, dec, radius, **kwargs):
        '''
        

        Parameters
        ----------
        ra : TYPE
            DESCRIPTION.
        dec : TYPE
            DESCRIPTION.
        radius : TYPE
            DESCRIPTION.
        **kwargs : TYPE
            DESCRIPTION.

        Returns
        -------
        result : TYPE
            DESCRIPTION.

        '''

        if 'data' in kwargs :
            data = kwargs['data']
        if 'conditions' in kwargs:
            conditions = kwargs['conditions']
    
        query = self._adqlWriter(ra, dec, radius, data, conditions)

        job = Gaia.launch_job_async(query)
        result = job.get_results()
        
        print("Sample number of sources: {:d}".format(len(query)))
        return result
    
    def getAstrometry(self, ra, dec, radius):
        '''
        

        Parameters
        ----------
        ra : TYPE
            DESCRIPTION.
        dec : TYPE
            DESCRIPTION.
        radius : TYPE
            DESCRIPTION.

        Returns
        -------
        astro_cluster : TYPE
            DESCRIPTION.

        '''        
        astrometry = 'source_id, ra, ra_error, dec, dec_error, parallax, parallax_error, pmra, pmra_error, pmdec, pmdec_error'
        
        query = self._adqlWriter(ra, dec, radius, data=astrometry)
        
        job = Gaia.launch_job_async(query)
        astro_cluster = job.get_results()
        print(" Sample number of sources: {:d}".format(len(astro_cluster)))
        
        return astro_cluster
    
    def getPhotometry(self, ra, dec, radius):
        '''
        

        Parameters
        ----------
        ra : TYPE
            DESCRIPTION.
        dec : TYPE
            DESCRIPTION.
        radius : TYPE
            DESCRIPTION.

        Returns
        -------
        photo_cluster : TYPE
            DESCRIPTION.

        '''
        photometry = 'source_id, bp_rp, phot_g_mean_mag, phot_bp_rp_excess_factor, teff_gspphot, ruwe, astrometric_excess_noise_sig'
        
        query = self._adqlWriter(ra, dec, radius, data=photometry)
        
        job = Gaia.launch_job_async(query)
        photo_cluster = job.get_results()
        
        print(" Sample number of sources: {:d}".format(len(photo_cluster)))
        
        return photo_cluster
    
    def getRV(self, ra, dec, radius, conditions=None):
        '''
        

        Parameters
        ----------
        ra : TYPE
            DESCRIPTION.
        dec : TYPE
            DESCRIPTION.
        radius : TYPE
            DESCRIPTION.
        conditions : TYPE, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        rv : TYPE
            DESCRIPTION.

        '''
        if isinstance(ra, u.Quantity) and isinstance(dec, u.Quantity):
            ra = str(ra/u.deg)
            dec= str(dec/u.deg)
        else:
            ra = str(ra)
            dec = str(dec)
        radius = str(radius)
        circle = ra + "," + dec + "," + radius
        query = "SELECT source_id, radial_velocity, radial_velocity_error \
                FROM gaiadr3.gaia_source \
                WHERE CONTAINS(POINT('ICRS',gaiadr3.gaia_source.ra,gaiadr3.gaia_source.dec),CIRCLE('ICRS'," + circle +"))=1"
                            
        job = Gaia.launch_job_async(query)
        rv = job.get_results()
        print(" Sample number of sources: {:d}".format(len(rv)))
        return rv
    
    def saveQuery(self, dat, name: str):
        '''
        

        Parameters
        ----------
        dat : TYPE
            DESCRIPTION.
        name : str
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        tn = self._tn()
        fold = self._checkPathExist(name.upper())
        path = fold+"/"+tn+'.txt'
        if isinstance(dat, Table)==False:
            data = Table(dat)
            
        data.write(path, format='ascii.tab')
        
        return (print(fold+tn+'.txt'))
    






        
