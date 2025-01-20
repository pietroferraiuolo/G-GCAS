"""
::module..ggcas._utility.sample

sample.py
=========


Author(s):
----------
- Pietro Ferraiuolo : Written in 2024

Description:
------------
This module provides a class for handling the query result sample. The class
unifies the Cluster object and the sample data obtained, in order to have a
compact object which has everything and from which quantities can be computed
easily.


"""

import pandas as _pd
from typing import List
from ggcas._cluster import Cluster
from typing import Optional, Union
from astropy.table import QTable, Table
from astropy import units as u

class Sample:
    """
    Class for handling the query result sample.

    It is an object with unifies the Cluster object and the sample data obtained,
    in order to have a compact object which has everythin and from which quantities
    can be computed easily.

    Parameters
    ----------
    gc : ggcas.cluster.Cluster
        Globular cluster object used for the query.
    sample : astropy.table.Table
        Table containing the retrieved sample's data.
    """
    def __init__(self, sample, gc:Optional[Union[Cluster,str]]=None):
        """The constructor"""
        self.gc     = gc if isinstance(gc, Cluster) else Cluster(gc)
        self.qinfo  = None
        self._sample = QTable.from_pandas(sample) if isinstance(sample, _pd.DataFrame) else sample
        self._table = None
        self._bckupSample = self._sample.copy()

    def __str__(self):
        """The string representation"""
        return self.gc.__str__()+'\n'+self._sample.__str__()

    def __repr__(self):
        """The representation"""
        return self.__get_repr()

    def __len__(self):
        """The length of the sample"""
        return len(self._sample['SOURCE_ID'])
    
    def __getitem__(self, key):
        """The item getter"""
        return self._sample[key]
    
    def __getattr__(self, attr):
        """The attribute getter"""
        if attr in self._sample.colnames:
            return self._sample[attr]
        else:
            raise AttributeError(f"'Sample' object has no attribute '{attr}'")
    
    def __setitem__(self, key, value):
        """The item setter"""
        self._sample[key] = value

    def __contains__(self, item):
        """The item checker"""
        return item in self._sample.colnames
    
    def __iter__(self):
        """The iterator"""
        return iter(self._sample.colnames)
    
    def __reversed__(self):
        """The reversed iterator"""
        return reversed(self._sample.colnames)
    
    @property
    def sample(self):
        """Returns the sample data"""
        return self._sample

    def computeDistance(self):
        """
        Computes the line-of-sight distance of each source in the sample from 
        Gaia, using the parallax.

        Results are stored in the 'rlos' and 'rlos_errors' columns.
        """
        from ggcas.functions import LosDistance
        los = LosDistance()
        data = [self._sample['parallax']]
        errs = [self._sample['parallax_error']]
        print('\nComputing Line-of-Sight Distance...\n')
        los.compute(data, errs)
        print('-'*30)
        self._sample['rlos'] = los.computed_values
        self._sample['rlos_errors'] = los.computed_errors
        return self._sample['rlos', 'rlos_errors'].info()
    
    def computeAngularSeparation(self):
        """
        Computes the angular separation of each source in the sample from the center
        of the cluster, using it's coordinates.

        Results are stored in the 'angsep' column.
        """
        from ggcas.functions import AngularSeparation
        angsep = AngularSeparation(self.gc.ra, self.gc.dec)
        data = [self._sample['ra'], self._sample['dec']]
        errs = [self._sample['ra_error'], self._sample['dec_error']]
        if 'ra_dec_corr' in self._sample.colnames:
            corr = [self._sample['ra_dec_corr']]
        else:
            corr = None
        print('\nComputing Angular Separation...\n')
        angsep.compute(data, errs, corr)
        print('-'*30)
        self._sample['angsep'] = angsep.computed_values
        self._sample['angsep_errors'] = angsep.computed_errors
        return self._sample[['angsep', 'angsep_errors']].info()
    
    def computeProjectedDistance(self):
        """
        Computes the projected distance of each source in the sample from the center
        of the cluster, using the angular separation and the cluster's distance from
        Earth.

        Results are stored in the 'r2d' and 'r2d_errors' columns.
        """
        from ggcas.functions import RadialDistance2D
        r2d = RadialDistance2D(self.gc.dist)
        if 'angsep' not in self._sample.colnames:
            self.computeAngularSeparation()
        data = [self._sample['angsep']]
        errs = [self._sample['angsep_errors']]
        print('\nComputing Projected Distance...\n')
        r2d.compute(data, errs)
        print('-'*30)
        self._sample['r2d'] = r2d.computed_values
        self._sample['r2d_errors'] = r2d.computed_errors
        return self._sample['r2d', 'r2d_errors'].info()
    
    def computeRadialDistance(self):
        """
        Computes the radial distance of each source in the sample from the
        center of the cluster, using it's coordinates.

        For this computation, which follows the formula:
        :math:`R_{3D} = \sqrt{R_{los}^2 + R_{2D}^2} = \sqrt{R_{los}^2 + (R_{gc} \\tan(\\theta_{x0}))^2}`

        where :math:`R_{los}` is the line-of-sight distance, :math:`R_{2D}` is 
        the projected distance from the center of the cluster, :math:`R_{gc}` is 
        the distance of the cluster from the Sun, and :math:`\\theta_{x0}` is 
        the angular separation of the source from the center of the cluster, if
        any of these quantities are not been computed yet, they will be computed first.

        In any case, the results are stored in the 'r3d' and 'r3d_errors' columns, as for the 
        additional needed quantities, if here computed, will be stored in the 'rlos'
        'rlos_errors', 'r2d', 'r2d_errors', 'angsep', 'angsep_errors' columns.
        """
        from ggcas.functions import RadialDistance3D
        if 'rlos' not in self._sample.colnames:
            self.computeDistance()
        if  'r2d' not in self._sample.colnames:
            self.computeProjectedDistance()
        data = [self._sample['rlos'], self._sample['r2d']]
        errs = [self._sample['rlos_errors'], self._sample['r2d_errors']]
        r3d = RadialDistance3D()
        print('\nComputing Radial Distance...\n')
        r3d.compute(data, errs)
        print('-'*30)
        self._sample['r3d'] = r3d.computed_values
        self._sample['r3d_errors'] = r3d.computed_errors
        return self._sample['r3d', 'r3d_errors'].info()

    def drop_columns(self, columns:list):
        """
        Drops the specified columns from the sample data.

        Parameters
        ----------
        columns : list
            List of column names to drop.
        """
        self._sample.remove_columns(columns)
    
    def info(self, *args):
        """Returns the info of the sample"""
        return self._sample.info(*args)

    def join(self, other, inplace:bool=False):
        """
        Joins the sample data with another sample data.

        Parameters
        ----------
        other : ggcas.Sample
            The other sample data to join with.
        inplace : bool
            If True, the operation is done in place, otherwise a new object is returned.

        Returns
        -------
        sample : ggcas.Sample
            The sample object containing the joined data.
        """
        sample = self._sample.to_pandas()
        other_sample = other.to_pandas()
        merged = sample.merge(other_sample, how='outer', indicator=True)
        merged_qtable = QTable.from_pandas(merged)
        if inplace:
            self._sample = merged_qtable
            self._merge_info = merged[['SOURCE_ID','_merge']]
            self.drop_columns(['_merge'])
            return merged_qtable
        else:
            new_sample = Sample(merged_qtable, self.gc)
            new_sample._merge_info = merged[['SOURCE_ID','_merge']]
            new_sample.drop_columns(['_merge'])
            return new_sample

    def to_pandas(self, overwrite:bool=False, *args, **kwargs):
        """
        Converts the sample (astropy.Table as default) to a pandas DataFrame.

        Parameters
        ----------
        *args : tuple
            Positional arguments to pass to the astropy.Table to_pandas method.
        **kwargs : dict
            Keyword arguments to pass to the astropy.Table to_pandas method.

        Returns
        -------
        df : pandas.DataFrame
            The DataFrame containing the sample data.
        """
        if isinstance(self._sample, (QTable, Table)):
            pandas_df = self._sample.to_pandas(*args, **kwargs)
            if overwrite:
                self._table = self._sample
                self._sample = pandas_df
                return self._sample.head(5)
            return pandas_df
        else:
            pass
    
    def to_table(self, *args):
        """
        Converts back the sample from a pandas.DataFrame into an astropy Table.

        Parameters
        ----------
        *args : tuple
            Positional arguments to pass to the astropy.Table constructor.

        Returns
        -------
        table : astropy.Table
            The table containing the sample data.
        """
        if not isinstance(self._sample, (QTable, Table)):
            if self._table is not None:
                if self._table.columns==self._sample.columns:
                    self._sample = self._table
                else:
                    self._sample = QTable.from_pandas(self._sample, *args)
            else:
                if isinstance(self._sample, _pd.DataFrame):
                    self._sample = QTable.from_pandas(self._sample, *args)
            return self._sample
        else: pass

    def to_numpy(self, columns:List[str]=None):
        """
        Converts the sample data to a numpy array.

        Returns
        -------
        arr : numpy.ndarray
            The numpy array containing the sample data.
        """
        if columns is not None:
            return self._sample[columns].to_pandas().to_numpy()
        else:
            return self._sample.to_pandas().to_numpy()
    
    def reset_sample(self):
        """Resets the sample to its original state"""
        self._sample = self._bckupSample.copy()

    def update_gc_params(self, **kwargs):
        """
        Updates the parameters of the cluster object.

        Parameters
        ----------
        **kwargs : dict
            The parameters to update.
        """
        for key in kwargs:
            if hasattr(self.gc, key):
                setattr(self.gc, key, kwargs[key])
            else:
                if not self.gc.id=='UntrackedData':
                    text = self.__get_repr()
                    text = text.split('\n')[5:]
                    ptxt = '\n'.join(text)
                else:
                    ptxt = ''
                raise AttributeError(f"'Cluster' object has no attribute '{key}'\n{ptxt}")
        return self.gc.__str__()


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
        i = 1
        for name in self._sample.colnames:
            stxt += name.lower()+' - '
            if i%5==0:
                stxt += '\n'
            i += 1
        stxt = stxt[:-3]
        return gctxt+stxt