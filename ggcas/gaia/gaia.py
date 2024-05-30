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
Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"
Gaia.ROW_LIMIT = -1


class GaiaQuery:
    
    def __init__(self):
        '''The Constructor'''
        self._path = "C:/Users/Er_da/Desktop/Poteff/data/query/"
        
    def _tn(self):
        '''
        

        Returns
        -------
        TYPE
            DESCRIPTION.

        '''
        return dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    
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
    
    def freeQuery(self, ra, dec, radius, data: list):
        '''
        

        Parameters
        ----------
        ra : TYPE
            DESCRIPTION.
        dec : TYPE
            DESCRIPTION.
        radius : TYPE
            DESCRIPTION.
        data : list
            DESCRIPTION.

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
        radius = str(radius)
        circle = ra + "," + dec + "," + radius
        select=""
        for i in range(len(data)-1):
            select += data[i]+', '
        select += data[len(data)-1]
        query = "SELECT "+ select +" \
                FROM gaiadr3.gaia_source \
                WHERE CONTAINS(POINT('ICRS',gaiadr3.gaia_source.ra,gaiadr3.gaia_source.dec),CIRCLE('ICRS'," + circle +"))=1"
        job = Gaia.launch_job_async(query)
        query = job.get_results()
        print("Sample number of sources: {:d}".format(len(query)))
        return query
    
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
        if isinstance(ra, u.Quantity) and isinstance(dec, u.Quantity):
            ra = str(ra/u.deg)
            dec= str(dec/u.deg)
        else:
            ra = str(ra)
            dec = str(dec)
        radius = str(radius)
        circle = ra + "," + dec + "," + radius
        query = "SELECT source_id, ra, ra_error, dec, dec_error, parallax, parallax_error, pmra, pmra_error, pmdec, pmdec_error \
                FROM gaiadr3.gaia_source \
                WHERE CONTAINS(POINT('ICRS',gaiadr3.gaia_source.ra,gaiadr3.gaia_source.dec),CIRCLE('ICRS'," + circle +"))=1"
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
        if isinstance(ra, u.Quantity) and isinstance(dec, u.Quantity):
            ra = str(ra/u.deg)
            dec= str(dec/u.deg)
        else:
            ra = str(ra)
            dec = str(dec)
        radius = str(radius)
        circle = ra + "," + dec + "," + radius
        query = "SELECT source_id, bp_rp, phot_g_mean_mag, phot_bp_rp_excess_factor, teff_gspphot \
                FROM gaiadr3.gaia_source \
                WHERE CONTAINS(POINT('ICRS',gaiadr3.gaia_source.ra,gaiadr3.gaia_source.dec),CIRCLE('ICRS'," + circle +"))=1"                     
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
    






        
