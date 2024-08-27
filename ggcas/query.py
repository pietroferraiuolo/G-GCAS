"""
Author(s)
---------
    - Pietro Ferraiuolo : Written in 2024

Description
-----------
Thi module contains the GaiaQuery class, which handles the ADQL language to make
queries for globular clusters data retrievement fast and easy.

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
from typing import Optional, Union
import os, configparser, numpy as np
from astropy.table import Table
from astropy import units as u
from astroquery.gaia import Gaia
from ggcas.cluster import Cluster
from ggcas.utility import folder_paths as fn, osutils as osu
from ggcas.utility.osutils import _timestamp
_QDATA = 'query_data.txt'
_QINFO = 'query_info.ini'

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

class Sample:
    """
    Class for better handling the query result sample.

    Parameters
    ----------
    gc : ggcas.cluster.Cluster
        Globular cluster object used for the query.
    sample : astropy.table.Table
        Table containing the retrieve sample's data.
    """
    def __init__(self, sample, gc:Cluster=None):
        """The constructor"""
        self.gc     = gc
        self.qinfo  = None
        self.sample = sample

    def __str__(self):
        """The string representation"""
        return self.gc.__str__()+'\n'+self.sample.__str__()

    def __repr__(self):
        """The representation"""
        return self.__get_repr()

    def __get_repr(self):
        """Gets the str representation"""
        if self.gc.id=='UntrackedData':
            gctxt=f"""Gaia data retrieved at coordinates
RA={self.gc.ra:.2f} DEC={self.gc.dec:.2f}
"""
        else:
            gctxt=f"""Data sample for cluster {self.gc.id}
"""
        stxt = '\nData Columns:\n'
        for name in self.sample.colnames:
            stxt += name.lower()+' - '
        stxt = stxt[:-3]
        return gctxt+stxt


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
    free_query:
        Perform an ADQL search into the Gaia catalogue with custom data to retrieve
        and conditions to apply.
    get_atrometry:
        A pre-constructed ADQL search into the Gaia catalogue with fixed data
        retrieved, which is the principal astrometric parameters, with the possibility
        to add personalized query conditions.
    get_photometry:
        A pre-constructed ADQL search into the Gaia catalogue with fixed data
        retrieved, which is the principal photometric parameters, with the possibility
        to add personalized query conditions.
    get_rv:
        A pre-constructed ADQL search into the Gaia catalogue with fixed data
        retrieved, which is the radial velociti parameter with its error, with
        the possibility to add personalized query conditions.

    How to Use it
    -------------
    Import the class and initialize it

    >>> from ggcas.query import GaiaQuery
    >>> dr3 = GaiaQuery()
    'Initialized with Gaia table: gaiadr3.gaia_source'

    To use a different Gaia catalogue simply initialize the class to it:

    >>> table = 'gaiadr2.gaia_source'
    >>> dr2 = GaiaQuery(gaia_table=table)
    'Initialized with Gaia table: gaiadr2.gaia_source'

    The queries, to work efficiently, require a 'ggcas.cluster.Cluster' object: so,
    let's take an example cluster, ngc104:

        >>> from ggcas.cluster import Cluster
        >>> gc = Cluster('ngc104')
        >>> gc
        <ggcas.cluster.Cluster object: NGC104>

    At this point, simply passing as argument the gc object to the query function:

        >>> data = dr3.get_astrometry(gc, radius=0.1) # Radius must be passed as deg
        INFO: Query finished. [astroquery.utils.tap.core]
        Sample number of sources: 45865
        >>> data
             SOURCE_ID              ra         ...        pmdec        pmdec_error
                                   deg         ...       mas / yr        mas / yr
               int64             float64       ...       float64         float32
        ------------------- ------------------ ... ------------------- -----------
        4689637445779054208   5.76563951251253 ...  -2.259232804525973  0.27563116
        4689638850232458368  5.874682871570303 ... -2.3177094407812033  0.17122078
                ...                ...         ...         ...             ...

    The methods also have the save option, to save both the result of the query
    and its information, such as the parameters used.

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
        self.last_result = None
        self.last_query = None
        print(f"Initialized with Gaia table: '{gaia_table}'")

    def __repr__(self):
        """The representation"""
        return self.__get_repr()

    def __str__(self):
        """The string representation"""
        return self.__get_str()

    def free_query(self, radius, gc:Cluster=None, save:bool=False, **kwargs):
        """
        This function allows to perform an ADQL search into the Gaia catalogue with
        personalized parameters, such as data to collect and conditions to apply.

        Parameters
        ----------
        radius : float
            Radius, in degrees, of the scan circle.
        gc : ggcas.cluster.Cluster
            Globular Cluster object created with the G-GCAS module.
        save : bool, optional
            Whether to save the obtained data with its information or not.
        **kwargs : additional optional arguments
            ra : float or str
                Right ascension coordinate for the centre of the scan (if no gc is provided).
            dec : float or str
                Declination coordinate for the centre of the scan (if no gc is provided)
            name : str
                String which provides the folder name where to save the data.
                Needed if no 'gc' object is supplied: if it is not given, and
                save was True, the data will be stored in the 'UntrackedData'
                folder.
            data : str or list of str
                List of parameters to retrieve, from the ones printed by ''.print_table()''.
                If this argument is missing, the only parameter retrieved is 'source_id'.
                Aliases:
                    'dat'
                    'data'
                    'params'
                    'parameters'
            conditions : str or list of str
                Listo of conditions on the parameters to apply upon scanning the
                archive. If no conditions are supplied, no conditions are applied.
                Aliases:
                    'cond'
                    'conds'
                    'condition'
                    'conditions'

        Returns
        -------
        job_result : astropy table
            Result of the async query, stored into an astropy table.

        """
        if gc is None:
            ra = kwargs.get('ra', None)
            dec = kwargs.get('dec', None)
            gc = Cluster(ra=ra, dec=dec)
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
            'Flag': {'Query': 'free'}
            }
        dat = osu.get_kwargs(('data', 'dat','params', 'parameters'), 'source_id', kwargs)
        self._queryInfo['Scan Info']['Data Acquired'],_ = self._formatCheck(dat, 'None')
        cond = osu.get_kwargs(('cond', 'conds', 'conditions', 'condition'), 'None', kwargs)
        if isinstance(cond, list):
            ccond = ''
            for c in range(len(cond)-1):
                ccond += cond[c]+', '
            ccond += cond[-1]
            self._queryInfo['Scan Info']['Conditions Applied'] = ccond
        else:
            self._queryInfo['Scan Info']['Conditions Applied'] = cond
        samp = self._run_query(savename, ra, dec, radius, dat, cond, save)
        sample = Sample(samp, gc=gc)
        sample.qinfo = self._queryInfo['Scan Info']
        return sample

    def get_astrometry(self, radius, gc:Cluster=None, save:bool=False, **kwargs):
        """
        A pre-constructed ADQL search into the Gaia catalogue with fixed data
        retrieved, which is the principal astrometric parameters, with the possibility
        to add personalized query conditions.

        The retrieved data is:
        'source_id, ra, ra_error, dec, dec_error, parallax, parallax_error, pmra,
        pmra_error, pmdec, pmdec_error'

        Parameters
        ----------
        radius : float
            Radius, in degrees, of the scan circle.
        gc : ggcas.cluster.Cluster
            Globular Cluster object created with the G-GCAS module.
        save : bool, optional
            Whether to save the obtained data with its information or not.
        **kwargs : additional optional arguments
            ra : float or str
                Right ascension coordinate for the centre of the scan (if no gc is provided).
            dec : float or str
                Declination coordinate for the centre of the scan (if no gc is provided)
            name : str
                String which provides the folder name where to save the data.
                Needed if no 'gc' object is supplied: if it is not given, and
                save was True, the data will be stored in the 'UntrackedData'
                folder.
            conditions : str or list of str
                Listo of conditions on the parameters to apply upon scanning the
                archive. If no conditions are supplied, no conditions are applied.
                Aliases:
                    'cond'
                    'conds'
                    'condition'
                    'conditions'

        Returns
        -------
        astro_cluster : astropy.Table
            Astropy table with  the query results.
        """
        if gc is None:
            ra = kwargs.get('ra', None)
            dec = kwargs.get('dec', None)
            gc = Cluster(ra=ra, dec=dec)
            savename = kwargs.get('name', 'UntrackedData')
        else:
            ra = gc.ra
            dec = gc.dec
            savename = gc.id
        astrometry = 'source_id, ra, ra_error, dec, dec_error, parallax, parallax_error, pmra, pmra_error, pmdec, pmdec_error'
        self._queryInfo = {
            'Scan Info': {
                'RA': ra,
                'DEC': dec,
                'Scan Radius': radius,
                'Data Acquired': astrometry
                },
            'Flag': {'Query': 'astrometry'}
            }
        cond = osu.get_kwargs(('cond', 'conds', 'conditions', 'condition'), 'None', kwargs)
        if isinstance(cond, list):
            ccond = ''
            for c in range(len(cond)-1):
                ccond += (cond[c]+', ')
            ccond += cond[-1]
            self._queryInfo['Scan Info']['Conditions Applied'] = ccond
        else:
            self._queryInfo['Scan Info']['Conditions Applied'] = cond
        astro_cluster = self._run_query(savename, ra, dec, radius, astrometry, cond, save)
        astro_sample = Sample(astro_cluster, gc=gc)
        astro_sample.qinfo = self._queryInfo['Scan Info']
        return astro_sample

    def get_photometry(self, radius, gc:Cluster=None, save:str=False, **kwargs):
        """
        A pre-constructed ADQL search into the Gaia catalogue with fixed data
        retrieved, which is the principal photometric parameters, with the possibility
        to add personalized query conditions.

        The retrieved data is:
        'source_id, bp_rp, phot_bp_mean_flux, phot_rp_mean_flux, phot_g_mean_mag,
        phot_bp_rp_excess_factor, teff_gspphot'

        Parameters
        ----------
        radius : float
            Radius, in degrees, of the scan circle.
        gc : ggcas.cluster.Cluster
            Globular Cluster object created with the G-GCAS module.
        save : bool, optional
            Whether to save the obtained data with its information or not.
        **kwargs : additional optional arguments
            ra : float or str
                Right ascension coordinate for the centre of the scan (if no gc is provided).
            dec : float or str
                Declination coordinate for the centre of the scan (if no gc is provided)
            name : str
                String which provides the folder name where to save the data.
                Needed if no 'gc' object is supplied: if it is not given, and
                save was True, the data will be stored in the 'UntrackedData'
                folder.
            conditions : str or list of str
                Listo of conditions on the parameters to apply upon scanning the
                archive. If no conditions are supplied, no conditions are applied.
                Aliases:
                    'cond'
                    'conds'
                    'condition'
                    'conditions'

        Returns
        -------
        photo_cluster : astropy.Table
            Astropy table with the results.
        """
        if gc is None:
            ra = kwargs.get('ra', None)
            dec = kwargs.get('dec', None)
            gc = Cluster(ra=ra, dec=dec)
            savename = kwargs.get('name', 'UntrackedData')
        else:
            ra = gc.ra
            dec = gc.dec
            savename = gc.id
        photometry = 'source_id, bp_rp, phot_bp_mean_flux, phot_rp_mean_flux, phot_g_mean_mag, phot_bp_rp_excess_factor, teff_gspphot'
        self._queryInfo = {
            'Scan Info': {
                'RA': ra,
                'DEC': dec,
                'Scan Radius': radius,
                'Data Acquired': photometry
                },
            'Flag': {'Query': 'photometry'}
            }
        cond = osu.get_kwargs(('cond', 'conds', 'conditions', 'condition'), 'None', kwargs)
        if isinstance(cond, list):
            ccond = ''
            for c in range(len(cond)-1):
                ccond += cond[c]+', '
            ccond += cond[-1]
            self._queryInfo['Scan Info']['Conditions Applied'] = ccond
        else:
            self._queryInfo['Scan Info']['Conditions Applied'] = cond
        photo_cluster = self._run_query(savename, ra, dec, radius, photometry, cond, save)
        phot_sample = Sample(photo_cluster, gc=gc)
        phot_sample.qinfo = self._queryInfo['Scan Info']
        return phot_sample

    def get_rv(self, radius, gc:Cluster=None, save:bool=False, **kwargs):
        """
        A pre-constructed ADQL search into the Gaia catalogue with fixed data
        retrieved, which is the radial velociti parameter with its error, with
        the possibility to add personalized query conditions.

        The retrieved data is:
        'source_id, radial_velocity, radial_velocity_error'

        Parameters
        ----------
        radius : float
            Radius, in degrees, of the scan circle.
        gc : ggcas.cluster.Cluster
            Globular Cluster object created with the G-GCAS module.
        save : bool, optional
            Whether to save the obtained data with its information or not.
        **kwargs : additional optional arguments
            ra : float or str
                Right ascension coordinate for the centre of the scan (if no gc is provided).
            dec : float or str
                Declination coordinate for the centre of the scan (if no gc is provided)
            name : str
                String which provides the folder name where to save the data.
                Needed if no 'gc' object is supplied: if it is not given, and
                save was True, the data will be stored in the 'UntrackedData'
                folder.
            conditions : str or list of str
                Listo of conditions on the parameters to apply upon scanning the
                archive. If no conditions are supplied, no conditions are applied.
                Aliases:
                    'cond'
                    'conds'
                    'condition'
                    'conditions'
        Returns
        -------
        rv_cluster : astropy.Table
            Astropy t able with te result.
        """
        if gc is None:
            ra = kwargs.get('ra', None)
            dec = kwargs.get('dec', None)
            gc = Cluster(ra=ra, dec=dec)
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
                'Data Acquired': rv
                },
            'Flag': {'Query': 'radvel'}
            }
        cond = osu.get_kwargs(('cond', 'conds', 'conditions', 'condition'), 'None', **kwargs)
        if isinstance(cond, list):
            ccond = ''
            for c in range(len(cond)-1):
                ccond += cond[c]+', '
            ccond += cond[-1]
            self._queryInfo['Scan Info']['Conditions Applied'] = ccond
        else:
            self._queryInfo['Scan Info']['Conditions Applied'] = cond
        rv_cluster = self._run_query(savename, ra, dec, radius, rv, cond, save)
        rv_sample = Sample(rv_cluster, gc=gc)
        rv_sample.qinfo = self._queryInfo['Scan Info']
        return rv_cluster

    def _run_query(self, gc_id, ra, dec, radius, data, cond, save):
        """
        The actual sub-routine which sends the query, checking if data already
        exists with the same query conditions, in which case loading it.

        Parameters
        ----------
        gc_id : str
            Folder name, used to check de existing data and to save the qued data
            (if save=True). Usually is the Cluster's name.'
        ra : str
            Right ascension of the centre of the scan.
        dec : str
            Declination of the centre of the scan.
        radius : str
            Radius of the circular scan.
        data : str
            The data to retrieve from the query.
        conditions : list of str
            The conditions to apply on the requested data.
        save : bool
            Wether to save or not the retrieved data.
        """
        check = self.__check_query_exists(gc_id)
        if check is False:
            query = self._adqlWriter(ra, dec, radius, data=data, conditions=cond)
            job = Gaia.launch_job_async(query)
            sample = job.get_results()
            print(f"Sample number of sources: {len(sample):d}")
            self.last_result = sample
            if save:
                self._saveQuery(sample, gc_id)
        else:
            print(f"""Found data with the same conditions for object {gc_id} in
{check[1]}.
Loading it...""")
            sample = osu.load_query(check[1])
            self.last_result = check[1]
            print(f"Sample number of sources: {len(sample):d}")
        return sample

    def _saveQuery(self, dat, name:str):
        """
        Routine for saving the query with its information, in the 'query_data.txt'
        and 'query_info.txt' files

        Parameters
        ----------
        dat : astropy.Table
            The astropy table containing all the retrieved data.
        name : str
            Where to save the data, usually the cluster's name.

        """
        config = configparser.ConfigParser()
        tn = _timestamp()
        fold = self._checkPathExist(name.upper())
        tnfold = os.path.join(fold, tn)
        os.mkdir(tnfold)
        data = os.path.join(tnfold, _QDATA)
        info = os.path.join(tnfold, _QINFO)
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
        Check if the path exists, and if not creates it.

        Parameters
        ----------
        dest : str
            the path to check.

        """
        self._fold = fn.CLUSTER_DATA_FOLDER(dest)
        if not os.path.exists(self._fold):
            os.makedirs(self._fold)
            print(f"Path '{self._fold}' did not exist. Created.")
        return self._fold

    def _formatCheck(self, data: Optional[Union[str,list]], conditions: Optional[Union[str,list]]):
        """
        Function to check and correct the format the 'data' and 'conditions'
        variables were imput with.

        Parameters
        ----------
        data : Optional[Union[str,list]]
            The data variable, containing all the requested parameters.
        conditions : Optional[Union[str,list]]
            The conditions variable, containing all the conditions to apply to the
            query.

        Returns
        -------
        dat : str
            The correct format for the data variable.
        cond : list of str
            The correct format for the conditions variable.

        """
        dat  = ''
        cond = ''
        if data is not None:
            if isinstance(data, list):
                for i in range(len(data)-1):
                    dat += data[i]+', '
                dat += data[-1]
            else: dat=data
        else: dat='source_id'
        if conditions != 'None':
            if isinstance(conditions, str):
                conditions = conditions.split(',')
                cond = 'AND '
                for i in range(len(conditions)-1):
                    cond += conditions[i]+"""
    AND """
                cond += conditions[-1]
            else:
                cond = 'AND '
                for i in range(len(conditions)-1):
                    cond += conditions[i]+"""
    AND """
                cond += conditions[-1]
        return dat, cond

    def _adqlWriter(self, ra, dec, radius, data, conditions):
        """
        This function writes the query, correctly formatting all the variables
        in order to be accepted by the GAIA ADQL language.

        Parameters
        ----------
        ra : str
            Right ascension.
        dec : str
            Declination.
        radius : str
            Scan radius.
        data : str
            Data to retrieve.
        conditions : list of str
            Conditions to apply.

        Returns
        -------
        query : str
            The full string to input the query with.

        """
        if isinstance(ra, u.Quantity):
            ra = ra / u.deg
        if isinstance(dec, u.Quantity):
            dec= dec / u.deg
        if isinstance(radius, u.Quantity):
            radius  = radius / u.deg
        circle  = f"{ra},{dec},{radius:.3f}"
        dat, cond = self._formatCheck(data, conditions)
        query = self._baseQ.format(data=dat, table=self._table, circle=circle, cond=cond)
        self.last_query = query
        return query

    def __check_query_exists(self, name):
        """
        Checks wether the requested query already exist saved for the Cluster.

        Parameters
        ----------
        name : str
            Folder where to search for saved data. Usually is the cluster's name

        Returns
        -------
        check : bool or tuple
            If no data was found to be compatible with the search, this is False.
            Either way it is a tuple, which first elemnt is True while the second
            is the complete file path to the corresponding saved data.
        """
        config = configparser.ConfigParser()
        try:
            tns = osu.tnlist(name)
        except FileNotFoundError:
            return False
        check = False
        for tn in tns:
            file_path = os.path.join(tn, _QINFO)
            if os.path.exists(file_path):
                config.read(file_path)
                try:
                    data_acquired = config['Scan Info']['Data Acquired']
                    conditions_applied = config['Scan Info']['Conditions Applied']
                    scan_radius = config['Scan Info']['Scan Radius']
                except KeyError as e:
                    print(f"Key error: {e}")
                    continue
                if (data_acquired == self._queryInfo['Scan Info']['Data Acquired'] and
                    conditions_applied == self._queryInfo['Scan Info']['Conditions Applied'] and
                    scan_radius == str(self._queryInfo['Scan Info']['Scan Radius'])):
                    check = (True, os.path.join(tn, _QDATA))
                    break
        return check

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
