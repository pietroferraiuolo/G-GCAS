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

    >>> from ggcas import query
    >>> gq = query.GaiaQuery()
    Initialized with Gaia table: 'gaiadr3.gaia_source'

To check all the available Gaia mission tables

    >>> query.available_tables()
    INFO: Retrieving tables... [astroquery.utils.tap.core]
    INFO: Parsing tables... [astroquery.utils.tap.core]
    INFO: Done. [astroquery.utils.tap.core]
    external.apassdr9
    external.catwise2020
    ...

To show information about the loaded table(s), one can get the table description
with the '__repr__' method, that is:

    >>> gq
    GAIADR3.GAIA_SOURCE
    -------------------
    This table has an entry for every Gaia observed source as published with this
    data release. It contains the basic source parameters, in their final state
    as processed by the Gaia Data Processing and Analysis Consortium from the raw
    data coming from the spacecraft. The table is complemented with others containing
    information specific to certain kinds of objects
    (e.g.~Solar--system objects, non--single stars, variables etc.) and value--added
    processing (e.g.~astrophysical parameters etc.). Further array data types
    (spectra, epoch measurements) are presented separately via Datalink resources.

While for a complete list of parameters within the table, print out the query
object:

    >>> print(gq)
    0 solution_id
    1 designation
    2 source_id
    3 random_index
    4 ref_epoch
    .
    .
    .
    150 ebpminrp_gspphot_upper
    151 libname_gspphot
"""
import os, configparser, numpy as np
from typing import Optional, Union
from astropy.table import Table
from astropy import units as u
from astroquery.gaia import Gaia
from ggcas.cluster import Cluster
from ggcas.utility import folder_paths as fn
from ggcas.utility.osutils import _timestamp
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

    def __repr__(self):
        """The representation"""
        return self.__get_repr()

    def __str__(self):
        """The string representation"""
        return self.__get_str()

    def free_query(self, gc:Cluster, radius, save:bool=False, **kwargs):
        """
        This method allows to perform an ADQL search into the Gaia catalogue with
        personalized parameters, such as data to collect and conditions to apply.

        Parameters
        ----------
        gc : ggcas.cluster.Cluster
            Globular Cluster object created with the G-GCAS module.
        radius : float | ArrayLike
            Radius, in degrees, of the scan circle.
        save : bool, optional
            Whether to save the file for the qued data or not.
        **kwargs : additional optional arguments
            ra: Right ascension coordinate for the centre of the scan (if no gc
                is provided).
            dec: Declination coordinate for the centre of the scan (if no gc is
                 provided)
            name: String which provides the folder name where to save the data.
                  Needed if no 'gc' object is supplied: if it is not given, and
                  save was True, the data will be stored in the 'UntrackedData'
                  folder.
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
        if gc is None:
            ra = kwargs.get('ra', None)
            dec = kwargs.get('dec', None)
            savename = kwargs.get('name', 'UntrackedData')
        else:
            ra = gc.ra
            dec = gc.dec
            savename = gc.id
        self._queryInfo = {
            'Scan Info': {
                'RA': ra,
                'DEC': dec,
                'Scan Radius': radius
                },
            'Quey Flag' : 'free'
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
        if save:
            self._saveQuery(result, savename)
        return result

    def get_astrometry(self, gc:Cluster, radius, save:bool=False, **kwargs):
        """


        Parameters
        ----------
        ra : float | ArrayLike
            Right ascension coordinate, in degrees, of the centre of the scan.
        dec : float | ArrayLike
            Declination coordinate, in degrees, of the centre of the scan..
        radius : float | ArrayLike
            Radius, in degrees, of the scan circle.
        save : bool, optional
            Whether to save the file or not. If the quesry is to be saved, the
            save argument must be the name of the globular cluster, so that the
            software knows the folder where to save, with a new tracking number.
        **kwargs : additional optional arguments
            ra: Right ascension coordinate for the centre of the scan (if no gc
                is provided).
            dec: Declination coordinate for the centre of the scan (if no gc is
                 provided)
            name: String which provides the folder name where to save the data.
                  Needed if no 'gc' object is supplied: if it is not given, and
                  save was True, the data will be stored in the 'UntrackedData'
                  folder.
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
        if gc is None:
            ra = kwargs.get('ra', None)
            dec = kwargs.get('dec', None)
            savename = kwargs.get('name', 'UntrackedData')
        else:
            ra = gc.ra
            dec = gc.dec
            savename = gc.id
        astrometry = 'source_id, ra, ra_error, dec, dec_error, parallax, \
            parallax_error, pmra, pmra_error, pmdec, pmdec_error'
        self._queryInfo = {
            'Scan Info': {
                'RA': ra,
                'DEC': dec,
                'Scan Radius': radius,
                'Data Acquired': astrometry.split(',')
                },
            'Quey Flag' : 'astrometry'
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
                self._saveQuery(astro_cluster, savename)
            else:
                raise TypeError(f"'save' was {save}, but must be a string. \
                                Specify the name of the object or of \
                                    the destination folder")
        return astro_cluster

    def get_photometry(self, gc:Cluster, radius, save:str=False, **kwargs):
        """


        Parameters
        ----------
        ra : float | ArrayLike
            Right ascension coordinate, in degrees, of the centre of the scan.
        dec : float | ArrayLike
            Declination coordinate, in degrees, of the centre of the scan..
        radius : float | ArrayLike
            Radius, in degrees, of the scan circle.
        save : bool, optional
            Whether to save the file or not. If the quesry is to be saved, the
            save argument must be the name of the globular cluster, so that the
            software knows the folder where to save, with a new tracking number.
        **kwargs : additional optional arguments
            ra: Right ascension coordinate for the centre of the scan (if no gc
                is provided).
            dec: Declination coordinate for the centre of the scan (if no gc is
                 provided)
            name: String which provides the folder name where to save the data.
                  Needed if no 'gc' object is supplied: if it is not given, and
                  save was True, the data will be stored in the 'UntrackedData'
                  folder.
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
        if gc is None:
            ra = kwargs.get('ra', None)
            dec = kwargs.get('dec', None)
            savename = kwargs.get('name', 'UntrackedData')
        else:
            ra = gc.ra
            dec = gc.dec
            savename = gc.id
        photometry = 'source_id, bp_rp, phot_bp_mean_flux, phot_rp_mean_flux, \
            phot_g_mean_mag, phot_bp_rp_excess_factor, teff_gspphot'
        self._queryInfo = {
            'Scan Info': {
                'RA': ra,
                'DEC': dec,
                'Scan Radius': radius,
                'Data Acquired': photometry.split(',')
                },
            'Quey Flag' : 'photometry'
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
                self._saveQuery(photo_cluster, savename)
            else:
                raise TypeError(f"'save' was {save}, but must be a string. \
                                Specify the name of the object or of \
                                    the destination folder")
        return photo_cluster

    def get_rv(self, gc:Cluster, radius, save:bool=False, **kwargs):
        """


        Parameters
        ----------
        ra : float | ArrayLike
            Right ascension coordinate, in degrees, of the centre of the scan.
        dec : float | ArrayLike
            Declination coordinate, in degrees, of the centre of the scan..
        radius : float | ArrayLike
            Radius, in degrees, of the scan circle.
        save : bool, optional
            Whether to save the file or not. If the quesry is to be saved, the
            save argument must be the name of the globular cluster, so that the
            software knows the folder where to save, with a new tracking number.
        **kwargs : additional optional arguments
            ra: Right ascension coordinate for the centre of the scan (if no gc
                is provided).
            dec: Declination coordinate for the centre of the scan (if no gc is
                 provided)
            name: String which provides the folder name where to save the data.
                  Needed if no 'gc' object is supplied: if it is not given, and
                  save was True, the data will be stored in the 'UntrackedData'
                  folder.
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
        if gc is None:
            ra = kwargs.get('ra', None)
            dec = kwargs.get('dec', None)
            savename = kwargs.get('name', 'UntrackedData')
        else:
            ra = gc.ra
            dec = gc.dec
            savename = gc.id
        rv = 'source_id, radial_velocity, radial_velocity_error'
        self._queryInfo = {
            'Scan Info': {
                'RA': ra,
                'DEC': dec,
                'Scan Radius': radius,
                'Data Acquired': rv.split(',')
                },
            'Quey Flag' : 'radvel'
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
                self._saveQuery(rv_cluster, savename)
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
        dat  = ''
        cond = ''
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

    def __load_table(self):
        """Load The table instanced table(s)"""
        if isinstance(self._table, list):
            table = np.zeros(len(self._table), dtype=object)
            for i, t in enumerate(self._table):
                table[i] = Gaia.load_table(t)
        else:
            table = Gaia.load_table(self._table)
        return table

    def __get_repr(self):
        """Get text for '__repr__' method"""
        table = self.__load_table()
        text = ''
        if isinstance(table, np.ndarray):
            for t in table:
                text += f"\n{t.name.upper()}\n"+'-'*len(t.name)+f"\n{t.description}\n"
        else :
            text = f"{table.name.upper()}\n"+'-'*len(table.name)+f"\n{table.description}"
        text += '\n \n<ggcas.query.GaiaQuery class>'
        return text

    def __get_str(self):
        """Get text for '__str__' method"""
        tables = self.__load_table()
        if isinstance(tables, np.ndarray):
            cols = np.zeros(tables.shape[0], dtype=list)
            text = ''
            for t in tables:
                text += f"{t.name.upper()}"+' '*10
            for i,table in enumerate(tables):
                cols[i] = []
                for columns in table.columns:
                    cols[i].append(f"{columns.name}")
            line=''
            for i in range(len(cols[0])):
                for t, n in enumerate(cols):
                    try:
                        tab = " "*(max([len(c) for c in (n)]) - len(n[i]) +2)
                        line += f"{n[i]}"+tab
                    except IndexError:
                        tab = " "*(len(tables[t].name)+10)
                        line += "--"+tab
                line += '\n'
            text += ('\n'+line)
        else:
            text = ""
            i=0
            for column in tables.columns:
                text += f"{i} {column.name}\n"
                i+=1
        return text
