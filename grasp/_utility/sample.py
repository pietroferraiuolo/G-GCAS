"""
::module..grasp._utility.sample

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
from astropy import units as _u
from typing import List as _List
from grasp._utility.cluster import Cluster as _Cluster
from typing import Optional as _Optional, Union as _Union
from astropy.table import QTable as _QTable, Table as _Table


class Sample:
    """
    Class for handling the query result sample.

    It is an object with unifies the Cluster object and the sample data obtained,
    in order to have a compact object which has everythin and from which quantities
    can be computed easily.

    Parameters
    ----------
    sample : astropy.table.Table
        Table containing the retrieved sample's data.
    gc : grasp.cluster.Cluster, optional
        Globular cluster object used for the query.
    """

    def __init__(self, sample, gc: _Optional[_Union[_Cluster, str]] = None):
        """The constructor"""
        self.gc = (gc if isinstance(gc, _Cluster) else _Cluster(gc)) if gc else None
        self.qinfo = None
        self._sample = (
            _QTable.from_pandas(sample) if isinstance(sample, _pd.DataFrame) else sample
        )
        self.__check_simulation()
        self._table = None
        self._bckupSample = self._sample.copy()

    def __str__(self):
        """The string representation"""
        if self._is_simulation:
            return f"Simulated data sample"+ "\n" + self._sample.__str__()
        else:
            return self.gc.__str__() + "\n" + self._sample.__str__()

    def __repr__(self):
        """The representation"""
        return self.__get_repr()

    def __len__(self):
        """The length of the sample"""
        return len(self._sample["SOURCE_ID"])

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


    def drop_columns(self, columns: list):
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
    

    def head(self, n: int = 5):
        """Returns the first n rows of the sample"""
        return self._sample.to_pandas().head(n)
    

    def describe(self):
        """Returns the description of the sample"""
        return self._sample.to_pandas().describe()


    def join(self, other, inplace: bool = False):
        """
        Joins the sample data with another sample data.

        Parameters
        ----------
        other : grasp.Sample
            The other sample data to join with.
        inplace : bool
            If True, the operation is done in place, otherwise a new object is returned.

        Returns
        -------
        sample : grasp.Sample
            The sample object containing the joined data.
        """
        sample = self._sample.to_pandas()
        other_sample = other.to_pandas()
        merged = sample.merge(other_sample, how="outer", indicator=True)
        merged_qtable = _QTable.from_pandas(merged)
        if inplace:
            self._sample = merged_qtable
            self._merge_info = merged[["SOURCE_ID", "_merge"]]
            self.drop_columns(["_merge"])
            return merged_qtable
        else:
            new_sample = Sample(merged_qtable, self.gc)
            new_sample._merge_info = merged[["SOURCE_ID", "_merge"]]
            new_sample.drop_columns(["_merge"])
            return new_sample


    def to_pandas(self, overwrite: bool = False, *args, **kwargs):
        """
        Converts the sample (`astropy.QTable` as default) to a pandas DataFrame.

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
        if isinstance(self._sample, (_QTable, _Table)):
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
        if not isinstance(self._sample, (_QTable, _Table)):
            if self._table is not None:
                if self._table.columns == self._sample.columns:
                    self._sample = self._table
                else:
                    self._sample = _QTable.from_pandas(self._sample, *args)
            else:
                if isinstance(self._sample, _pd.DataFrame):
                    self._sample = _QTable.from_pandas(self._sample, *args)
            return self._sample
        else:
            pass


    def to_numpy(self, columns: _List[str] = None):
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
        if self._is_simulation:
            return "This is a simulation data sample. No GC available."
        for key in kwargs:
            if hasattr(self.gc, key):
                setattr(self.gc, key, kwargs[key])
            else:
                if not self.gc.id == "UntrackedData":
                    text = self.__get_repr()
                    text = text.split("\n")[5:]
                    ptxt = "\n".join(text)
                else:
                    ptxt = ""
                raise AttributeError(
                    f"'Cluster' object has no attribute '{key}'\n{ptxt}"
                )
        return self.gc.__str__()
    

    def apply_conditions(self, conditions: dict, inplace: bool = False):
        """
        Applies conditions to the sample data.

        Parameters
        ----------
        conditions : dict
            The conditions to apply.
        inplace : bool
            If True, the operation is done in place, otherwise a new object is returned.

        Returns
        -------
        sample : grasp.Sample
            The sample object containing the filtered data.

        How to Use
        ----------
        The correct use of the method resides completely on how the conditions dictionary
        is structured. The dictionary must have the column name as the key and the condition
        including the logical operation:

        >>> conditions = {
        ...     "parallax": "> -5",
        ...     "parallax": "< 5"
        ... }
        >>> sample.apply_conditions(conditions)
        """
        sample = self._sample.copy()
        conds = []
        for k,v in conditions.items():
            conds.append(f"(sample.{k} {v})")
        conds = " & ".join(conds)
        sample = sample[eval(conds)]
        if inplace:
            self._sample = sample
            return self.__repr__()
        else:
            return Sample(sample, self.gc)


    def __check_simulation(self):
        """Check wether the data the sample has been instanced with is
        a simulation or real data"""
        sim_a = [
            'Mass_[Msun]',
            'x_[pc]',
            'y_[pc]',
            'z_[pc]',
            'vx_[km/s]',
            'vy_[km/s]',
            'vz_[km/s]'
        ]
        if all(a == b for a, b in zip(sim_a, self.__iter__())):
            self.qinfo = "McLuster Simulation"
            self._sample['M'] = self._sample['Mass_[Msun]'] * _u.Msun
            self._sample["x"] = self._sample["x_[pc]"] * _u.pc
            self._sample["y"] = self._sample["y_[pc]"] * _u.pc
            self._sample["z"] = self._sample["z_[pc]"] * _u.pc
            self._sample["vx"] = self._sample["vx_[km/s]"] * _u.km / _u.s
            self._sample["vy"] = self._sample["vy_[km/s]"] * _u.km / _u.s
            self._sample["vz"] = self._sample["vz_[km/s]"] * _u.km / _u.s
            self.drop_columns(sim_a)
            self._is_simulation = True
        else:
            self._is_simulation = False
            return


    def __get_repr(self):
        """Gets the str representation"""
        if self._is_simulation:
            gctxt = f"""Simulated data sample"""
        elif self.gc.id == "UntrackedData":
            gctxt = f"""Gaia data retrieved at coordinates
RA={self.gc.ra:.2f} DEC={self.gc.dec:.2f}
"""
        else:
            gctxt = f"""Data sample for cluster {self.gc.id}
"""
        stxt = "\nData Columns:\n"
        i = 1
        for name in self._sample.colnames:
            stxt += name.lower() + " - "
            if i % 5 == 0:
                stxt += "\n"
            i += 1
        stxt = stxt[:-3]
        return gctxt + stxt
