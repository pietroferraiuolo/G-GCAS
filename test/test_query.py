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
import configparser

class GaiaQuery:
    
    def __init__(self, gaia_table: Optional[Union[str, list]] = "gaiadr3.gaia_source"):
        '''The Constructor
        
        Parameters
        ----------
        gaia_table: str | list
            DESCRIPTION.
        
        '''
        Gaia.MAIN_GAIA_TABLE = gaia_table
        Gaia.ROW_LIMIT = -1
        
        self._table     = gaia_table
        self._path      = os.environ['PYGCASCONF']
        self._fold      = None
        self._queryInfo = {}
        
        self._baseQ     = """SELECT {data}
                FROM {table} 
                WHERE CONTAINS(POINT('ICRS',gaiadr3.gaia_source.ra,gaiadr3.gaia_source.dec),CIRCLE('ICRS',{circle}))=1
                  {cond}
                """
        self._joinQ     = """SELECT {data}
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

        self._fold = os.path.join(self._path+'\query', dest.upper())
        if not os.path.exists(self._fold):
            os.makedirs(self._fold)
            print(f"Path '{self._fold}' did not exist. Created.")
            return self._fold
        else:
            return self._fold
        
    def _formatCheck(self, data: Optional[Union[str,list]], conditions: Optional[Union[str,list]]):
        '''
        

        Parameters
        ----------
        data : Optional[Union[str,list]]
            DESCRIPTION.
        conditions : Optional[Union[str,list]]
            DESCRIPTION.

        Returns
        -------
        dat : TYPE
            DESCRIPTION.
        cond : TYPE
            DESCRIPTION.

        '''
        dat     = ''
        cond    = ''
        if data is not None:
            if isinstance(data, list):
                for i in range(len(data)-1):
                    dat += data[i]+', '
                dat += data[len(data)-1]
            else: dat=data
        else: dat='source_id'
        
        if conditions is not None:
            if isinstance(conditions, str):
                conditions = conditions.split(',')
            cond = '  AND '
            for i in range(len(conditions)-1):
                cond += conditions[i]+"""
            AND """
            cond += conditions[len(conditions)-1]
        
        return dat, cond
        
    def _adqlWriter(self, ra, dec, radius, data, conditions):
        '''
        

        Parameters
        ----------
        ra : TYPE
            DESCRIPTION.
        dec : TYPE
            DESCRIPTION.
        radius : TYPE
            DESCRIPTION.
        data : TYPE
            DESCRIPTION.
        conditions : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        if isinstance(ra, u.Quantity) and isinstance(dec, u.Quantity):
            ra = str(ra/u.deg)
            dec= str(dec/u.deg)
        else:
            ra = str(ra)
            dec = str(dec)
            
        if isinstance(radius, u.Quantity):
            radius  = str(radius/u.deg)
        else:
            radius = str(radius)
        
        circle  = ra + "," + dec + "," + radius

        dat, cond = self._formatCheck(data, conditions)
        
        query = self._baseQ.format(data=dat, table=self._table, circle=circle, cond=cond)
        
        return query
        
    def printTable(self, dump = False):
        '''
        

        Parameters
        ----------
        dump : TYPE, optional
            DESCRIPTION. The default is False.

        Returns
        -------
        TYPE
            DESCRIPTION.

        '''
        table = Gaia.load_table(self._table)
        print(table.description)
        print('')
        i=0
        for columns in table.columns:
            print(i, columns.name)
            i+=1
            
        if dump==True:
            return table.columns
        
    
    def freeQuery(self, ra, dec, radius, save: str = False, **kwargs):
        '''
        

        Parameters
        ----------
        ra : TYPE
            DESCRIPTION.
        dec : TYPE
            DESCRIPTION.
        radius : TYPE
            DESCRIPTION.
        save : str, optional
            DESCRIPTION. The default is False.
        **kwargs : TYPE
            DESCRIPTION.

        Raises
        ------
        TypeError
            DESCRIPTION.

        Returns
        -------
        result : TYPE
            DESCRIPTION.

        '''
        self._queryInfo = {
            'Scan Info': {
                'RA': ra,
                'DEC': dec,
                'Scan Radius': radius
                }
            }
        
        if 'data' in kwargs :
            dat = kwargs['data']
            dat, _= self._formatCheck(dat, None) #gives a string
            self._queryInfo['Scan Info']['Data Acquired'] = dat.split(',')
        else: 
            dat='source_id'
            self._queryInfo['Scan Info']['Data Acquired'] = dat
            
        if 'conditions' in kwargs:
            cond = kwargs['conditions']
            if isinstance(cond, str):
                self._queryInfo['Scan Info']['Conditions Applied'] = cond.split(',')
            else:
                self._queryInfo['Scan Info']['Conditions Applied'] = cond
        else: 
            cond=None
            self._queryInfo['Scan Info']['Conditions Applied'] = 'None'
    
        query = self._adqlWriter(ra, dec, radius, dat, cond)

        job = Gaia.launch_job_async(query)
        result = job.get_results()
        
        print("Sample number of sources: {:d}".format(len(result)))
        
        if save is not False:
            if isinstance(save, str):
                self._saveQuery(result, save)
            else:
                raise TypeError("'save' was {save}, but must be a string. Please specify the name of the object or of the destination folder".format(save=type(save)))
        
        return result
    
    def getAstrometry(self, ra, dec, radius, save: str = False, **kwargs):
        '''
        

        Parameters
        ----------
        ra : TYPE
            DESCRIPTION.
        dec : TYPE
            DESCRIPTION.
        radius : TYPE
            DESCRIPTION.
        save : str, optional
            DESCRIPTION. The default is False.
        **kwargs : TYPE
            DESCRIPTION.

        Raises
        ------
        TypeError
            DESCRIPTION.

        Returns
        -------
        astro_cluster : TYPE
            DESCRIPTION.

        '''    
        astrometry = 'source_id, ra, ra_error, dec, dec_error, parallax, parallax_error, pmra, pmra_error, pmdec, pmdec_error'
        
        self._queryInfo = {
            'Scan Info': {
                'RA': ra,
                'DEC': dec,
                'Scan Radius': radius,
                'Data Acquired': astrometry.split(',')
                }
            }
        
        if 'conditions' in kwargs:
            cond = kwargs['conditions']
            if isinstance(cond, str):
                self._queryInfo['Scan Info']['Conditions Applied'] = cond.split(',')
            else:
                self._queryInfo['Scan Info']['Conditions Applied'] = cond
        else: 
            cond=None
            self._queryInfo['Scan Info']['Conditions Applied'] = 'None'
                    
        query = self._adqlWriter(ra, dec, radius, data=astrometry, conditions=cond)
        
        job = Gaia.launch_job_async(query)
        astro_cluster = job.get_results()
        print(" Sample number of sources: {:d}".format(len(astro_cluster)))
        
        if save is not False:
            if isinstance(save, str):
                self._saveQuery(astro_cluster, save)
            else:
                raise TypeError("'save' was {save}, but must be a string. Please specify the name of the object or of the destination folder".format(save=type(save)))
        
        return astro_cluster
    
    def getPhotometry(self, ra, dec, radius, save: str = False, **kwargs):
        '''
        

        Parameters
        ----------
        ra : TYPE
            DESCRIPTION.
        dec : TYPE
            DESCRIPTION.
        radius : TYPE
            DESCRIPTION.
        save : str, optional
            DESCRIPTION. The default is False.
        **kwargs : TYPE
            DESCRIPTION.

        Raises
        ------
        TypeError
            DESCRIPTION.

        Returns
        -------
        photo_cluster : TYPE
            DESCRIPTION.

        '''
        photometry = 'source_id, bp_rp, phot_bp_mean_flux, phot_rp_mean_flux, phot_g_mean_mag, phot_bp_rp_excess_factor, teff_gspphot'

        self._queryInfo = {
            'Scan Info': {
                'RA': ra,
                'DEC': dec,
                'Scan Radius': radius,
                'Data Acquired': photometry.split(',')
                }
            }
        
        if 'conditions' in kwargs:
            cond = kwargs['conditions']
            if isinstance(cond, str):
                self._queryInfo['Scan Info']['Conditions Applied'] = cond.split(',')
            else:
                self._queryInfo['Scan Info']['Conditions Applied'] = cond
        else: 
            cond=None
            self._queryInfo['Scan Info']['Conditions Applied'] = cond
            
        query = self._adqlWriter(ra, dec, radius, data=photometry, conditions=cond)    
        
        job = Gaia.launch_job_async(query)
        photo_cluster = job.get_results()
        
        print(" Sample number of sources: {:d}".format(len(photo_cluster)))
        
        if save is not False:
            if isinstance(save, str):
                self._saveQuery(photo_cluster, save)
            else:
                raise TypeError("'save' was {save}, but must be a string. Please specify the name of the object or of the destination folder".format(save=type(save)))
        
        return photo_cluster
    
    def getRV(self, ra, dec, radius, save: str = False, **kwargs):
        '''
        

        Parameters
        ----------
        ra : TYPE
            DESCRIPTION.
        dec : TYPE
            DESCRIPTION.
        radius : TYPE
            DESCRIPTION.
        save : str, optional
            DESCRIPTION. The default is False.
        **kwargs : TYPE
            DESCRIPTION.

        Raises
        ------
        TypeError
            DESCRIPTION.

        Returns
        -------
        rv_cluster : TYPE
            DESCRIPTION.

        '''
        rv = 'source_id, radial_velocity, radial_velocity_error'
        
        self._queryInfo = {
            'Scan Info': {
                'RA': ra,
                'DEC': dec,
                'Scan Radius': radius,
                'Data Acquired': rv.split(',')
                }
            }
        
        if 'conditions' in kwargs:
            cond = kwargs['conditions']
            if isinstance(cond, str):
                self._queryInfo['Scan Info']['Conditions Applied'] = cond.split(',')
            else:
                self._queryInfo['Scan Info']['Conditions Applied'] = cond
        else: 
            cond=None
            self._queryInfo['Scan Info']['Conditions Applied'] = cond
        
        query = self._adqlWriter(ra, dec, radius, data=rv, conditions=cond)
                            
        job = Gaia.launch_job_async(query)
        rv_cluster = job.get_results()
        print(" Sample number of sources: {:d}".format(len(rv_cluster)))
        
        if save is not False:
            if isinstance(save, str):
                self._saveQuery(rv_cluster, save)
            else:
                raise TypeError("'save' was {save}, but must be a string. Please specify the name of the object or of the destination folder".format(save=type(save)))
                
        return rv_cluster
    
    def _saveQuery(self, dat, name: str):
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
        config = configparser.ConfigParser()
        
        tn = self._tn()
        fold = self._checkPathExist(name.upper())
        
        path = os.path.join(fold, (tn+'.txt'))
        info = os.path.join(fold, (tn+'.ini'))
        
        if isinstance(dat, Table)==False:
            dat = Table(dat)
            
        dat.write(path, format='ascii.tab')
        
        for section, options in self._queryInfo.items():
            config[section] = options
        
        with open(info, 'w') as configfile:
            config.write(configfile)
        
        print(path)
        print(info)
