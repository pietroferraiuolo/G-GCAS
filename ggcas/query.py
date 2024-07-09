"""
Author(s)
---------
    - Pietro Ferraiuolo : Written in 2024

Description
-----------
Thi module contains the GaiaQuery class, which handles the ADQL language to mak
e queries for globular clusters data retrievement fast and easy. 

How to Use it
-------------
After importing the module, initialize the class with a table (default is the 
Gaia data release 3 table)
    
    >>> from ggcas.query import GaiaQuery
    >>> gq = GaiaQuery()
    Initialized with Gaia table: 'gaiadr3.gaia_source'
    
To check all the available Gaia mission tables

    >>> gq.available_tables()
    INFO: Retrieving tables... [astroquery.utils.tap.core]
    INFO: Parsing tables... [astroquery.utils.tap.core]
    INFO: Done. [astroquery.utils.tap.core]
    external.apassdr9
    external.catwise2020
    ...
"""
import os
import configparser
from typing import Optional, Union
from astropy.table import Table
from astropy import units as u
from astroquery.gaia import Gaia
from ggcas.utility import folder_paths as fn
from ggcas.utility.utils import _timestamp
QDATA = 'query_data.txt'
QINFO = 'query_info.ini'

def available_tables(key:str=None):
    """
    Prints out the complete list of data tables present in the Gaia archive.

    Parameters
    ----------
    key : str, optional
        A key used to restrict the printed tables. As example, if 
        >>> key = 'gaiadr3'
        then only tables relative to the complete 3th data release will be printed
        out. Default is None, meaning all the tables will be printed.
    """
    tables = Gaia.load_tables(only_names=True)
    if key is not None:
        for table in tables:
            name = table.name
            if key in name:
                print(name)
    else:
        for table in tables:
            print(table.name)

class GaiaQuery:
    """
    Classs for the Gaia Query language execution. 
    
    Description
    -----------
    With this class, it is possible to easily perform async queries and retriev
    e data from the ESA/GAIA catalogue. It is possible to use different data re
    leases by loading different data tables in the initialization of the class.
    
    Methods
    -------
    print_table:
        Print the loaded data table information.
    free_query:
        Perform an ADQL search into the Gaia catalogue with full personalized 
        parameters.
    get_atrometry:
        D
    get_photometry:
        D
    get_rv:
        D

    How to Use it
    -------------
    Import the class and initialize it
    
    >>> from ggcas.query import GaiaQuery
    >>> gq = GaiaQuery()
    'Initialized with Gaia table: gaiadr3.gaia_source'
    
    To use a different Gaia catalogue simply initialize the class to it:
        
    >>> table = 'gaiadr2.gaia_source'
    >>> gq = GaiaQuery(gaia_table=table)
    'Initialized with Gaia table: gaiadr2.gaia_source'
    """
    def __init__(self, gaia_table: Optional[Union[str, list]] = "gaiadr3.gaia_source"):
        """
        The Constructor

        Parameters
        ----------
        gaia_table : str or list of str, optional
            Gaia table(s) to initialize the class with. The default is the 3th 
            Gaia data release "gaiadr3.gaia_source".
        """
        Gaia.MAIN_GAIA_TABLE = gaia_table
        Gaia.ROW_LIMIT = -1
        self._table     = gaia_table
        self._path      = fn.BASE_DATA_PATH
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
        print(f"Initialized with Gaia table: '{gaia_table}'")

    def print_table(self, dump:bool = False):
        """
        Print the loaded data table information. If dump is false, each column
        name is printed too, with its index number.

        Parameters
        ----------
        dump : bool, optional
            Option for storing the table's cloumns into a variable. The default
            is False, and the columns will be printed with a numbered index.

        Returns
        -------
        table : list of astroquery TapColumn
            A list in which each element is a column of the loaded Gaia data table.
            Printing and entry will display all its information.
        """
        table = Gaia.load_table(self._table)
        print(table.description)
        print('')
        i=0
        for columns in table.columns:
            print(i, columns.name)
            i+=1
        if dump:
            return table.columns

    def free_query(self, ra, dec, radius, save: str = False, **kwargs):
        """
        This method allows to perform an ADQL search into the Gaia catalogue with
        personalized parameters, such as data to collect and conditions to apply.

        Parameters
        ----------
        ra : float | ArrayLike
            Right ascension coordinate, in degrees, of the centre of the scan.
        dec : float | ArrayLike
            Declination coordinate, in degrees, of the centre of the scan..
        radius : float | ArrayLike
            Radius, in degrees, of the scan circle.
        save : str, optional
            Whether to save the file or not. If the quesry is to be saved, the 
            save argument must be the name of the globular cluster, so that the
            software knows the folder where to save, with a new tracking number.
        **kwargs : additional optional arguments
            data: str or list of str
                List of parameters to retrieve, from the ones printed by ''.print_table()''.
                If this argument is missing, the only parameter retrieved is 'source_id'.
            conditions: str or list of str
                Listo of conditions on the parameters to apply upon scanning the
                archive. If no conditions are supplied, no conditions are applied.

        Returns
        -------
        result : astropy table
            Result of the async query, stored into an astropy table.

        Raises
        ------
        TypeError
            Raised if save is not a string.

        Example
        -------
        """
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
        print(f"Sample number of sources: {(len(result)):d}")
        if save is not False:
            if isinstance(save, str):
                self._saveQuery(result, save)
            else:
                raise TypeError(f"'save' was {save}, but must be a string. \
                                Specify the name of the object or of the\
                                    destination folder")
        return result

    def get_astrometry(self, ra, dec, radius, save:str=False, **kwargs):
        """
        

        Parameters
        ----------
        ra : float | ArrayLike
            Right ascension coordinate, in degrees, of the centre of the scan.
        dec : float | ArrayLike
            Declination coordinate, in degrees, of the centre of the scan..
        radius : float | ArrayLike
            Radius, in degrees, of the scan circle.
        save : str, optional
            Whether to save the file or not. If the quesry is to be saved, the 
            save argument must be the name of the globular cluster, so that the
            software knows the folder where to save, with a new tracking number.
        **kwargs : additional optional arguments
            conditions: str or list of str
                Listo of conditions on the parameters to apply upon scanning the
                archive. If no conditions are supplied, no conditions are applied.

        Returns
        -------
        astro_cluster : TYPE
            DESCRIPTION.

        Raises
        ------
        TypeError
            Raised if save is not a string.

        """
        astrometry = 'source_id, ra, ra_error, dec, dec_error, parallax, \
            parallax_error, pmra, pmra_error, pmdec, pmdec_error'
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
        print(f"Sample number of sources: {len(astro_cluster):d}")
        if save is not False:
            if isinstance(save, str):
                self._saveQuery(astro_cluster, save)
            else:
                raise TypeError(f"'save' was {save}, but must be a string. \
                                Specify the name of the object or of \
                                    the destination folder")
        return astro_cluster

    def get_photometry(self, ra, dec, radius, save: str = True, **kwargs):
        """
        

        Parameters
        ----------
        ra : float | ArrayLike
            Right ascension coordinate, in degrees, of the centre of the scan.
        dec : float | ArrayLike
            Declination coordinate, in degrees, of the centre of the scan..
        radius : float | ArrayLike
            Radius, in degrees, of the scan circle.
        save : str, optional
            Whether to save the file or not. If the quesry is to be saved, the 
            save argument must be the name of the globular cluster, so that the
            software knows the folder where to save, with a new tracking number.
        **kwargs : additional optional arguments
            conditions: str or list of str
                Listo of conditions on the parameters to apply upon scanning the
                archive. If no conditions are supplied, no conditions are applied.

        Returns
        -------
        photo_cluster : TYPE
            DESCRIPTION.

        Raises
        ------
        TypeError
            Raised if save is not a string.
        """
        photometry = 'source_id, bp_rp, phot_bp_mean_flux, phot_rp_mean_flux, \
            phot_g_mean_mag, phot_bp_rp_excess_factor, teff_gspphot'
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
        print(f"Sample number of sources: {len(photo_cluster):d}")
        if save is not False:
            if isinstance(save, str):
                self._saveQuery(photo_cluster, save)
            else:
                raise TypeError(f"'save' was {save}, but must be a string. \
                                Specify the name of the object or of \
                                    the destination folder")
        return photo_cluster

    def get_rv(self, ra, dec, radius, save: str = True, **kwargs):
        """
        

        Parameters
        ----------
        ra : float | ArrayLike
            Right ascension coordinate, in degrees, of the centre of the scan.
        dec : float | ArrayLike
            Declination coordinate, in degrees, of the centre of the scan..
        radius : float | ArrayLike
            Radius, in degrees, of the scan circle.
        save : str, optional
            Whether to save the file or not. If the quesry is to be saved, the 
            save argument must be the name of the globular cluster, so that the
            software knows the folder where to save, with a new tracking number.
        **kwargs : additional optional arguments
            conditions: str or list of str
                Listo of conditions on the parameters to apply upon scanning the
                archive. If no conditions are supplied, no conditions are applied.
        Returns
        -------
        rv_cluster : TYPE
            DESCRIPTION.

        Raises
        ------
        TypeError
            Raised if save is not a string.
        """
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
        print(f"Sample number of sources: {len(rv_cluster):d}")
        if save is not False:
            if isinstance(save, str):
                self._saveQuery(rv_cluster, save)
            else:
                raise TypeError(f"'save' was {save}, but must be a string. \
                                Specify the name of the object or of the destination folder")
        return rv_cluster

    def _saveQuery(self, dat, name: str):
        """
        

        Parameters
        ----------
        dat : TYPE
            DESCRIPTION.
        name : str
            DESCRIPTION.

        """
        config = configparser.ConfigParser()
        tn = _timestamp()
        fold = self._checkPathExist(name.upper())
        tnfold = os.path.join(fold, tn)
        os.mkdir(tnfold)
        data = os.path.join(tnfold, QDATA)
        info = os.path.join(tnfold, QINFO)
        if isinstance(dat, Table) is False:
            dat = Table(dat)
        dat.write(data, format='ascii.tab')
        for section, options in self._queryInfo.items():
            config[section] = options
        with open(info, 'w', encoding='UTF-8') as configfile:
            config.write(configfile)
        print(data)
        print(info)

    def _checkPathExist(self, dest: str):
        """
        

        Parameters
        ----------
        dest : str
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        self._fold = fn.CLUSTER_DATA_FOLDER(dest)
        if not os.path.exists(self._fold):
            os.makedirs(self._fold)
            print(f"Path '{self._fold}' did not exist. Created.")
        return self._fold

    def _formatCheck(self, data: Optional[Union[str,list]], conditions: Optional[Union[str,list]]):
        """
        

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

        """
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
        """
        

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
        query : TYPE
            DESCRIPTION.

        """
        if isinstance(ra, u.Quantity) and isinstance(dec, u.Quantity):
            ra = str(ra / u.deg)
            dec= str(dec / u.deg)
        else:
            ra = str(ra)
            dec = str(dec)
        if isinstance(radius, u.Quantity):
            radius  = str(radius / u.deg)
        else:
            radius = str(radius)
        circle  = ra + "," + dec + "," + radius
        dat, cond = self._formatCheck(data, conditions)
        query = self._baseQ.format(data=dat, table=self._table, circle=circle, cond=cond)
        return query

    def __check_query_exists(self):
        
        return
    