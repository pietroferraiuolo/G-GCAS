from typing import Optional, Union, List, Tuple, Dict
import pandas
from astropy.table import Table
from ggcas._cluster import Cluster

def available_tables(key: Optional[str] = None) -> None: ...

class GaiaQuery:
    def __init__(
        self, gaia_table: Optional[Union[str, List[str]]] = "gaiadr3.gaia_source"
    ) -> None: ...

    def __repr__(self) -> str: ...

    def __str__(self) -> str: ...

    def free_query(
        self,
        radius: float,
        gc: Optional[Union[Cluster, str]] = None,
        save: bool = False,
        **kwargs
    ) -> "_Sample": ...

    def get_astrometry(
        self,
        radius: float,
        gc: Optional[Union[Cluster, str]] = None,
        save: bool = False,
        **kwargs
    ) -> "_Sample": ...

    def get_photometry(
        self,
        radius: float,
        gc: Optional[Union[Cluster, str]] = None,
        save: bool = False,
        **kwargs
    ) -> "_Sample": ...

    def get_rv(
        self,
        radius: float,
        gc: Optional[Union[Cluster, str]] = None,
        save: bool = False,
        **kwargs
    ) -> "_Sample": ...

    def _run_query(
        self,
        gc_id: str,
        ra: float,
        dec: float,
        radius: float,
        data: Union[str, List[str]],
        cond: Union[str, List[str]],
        save: bool,
    ) -> Table: ...

    def _saveQuery(self, dat: Table, name: str) -> None: ...

    def _checkPathExist(self, dest: str) -> None: ...

    def _formatCheck(
        self,
        data: Optional[Union[str, List[str]]],
        conditions: Optional[Union[str, List[str]]],
    ) -> Tuple[str, str]: ...

    def _adqlWriter(
        self,
        ra: float,
        dec: float,
        radius: float,
        data: Union[str, List[str]],
        conditions: Union[str, List[str]],
    ) -> str: ...

    def __check_query_exists(self, name: str) -> bool: ...

    def __load_table(self) -> Table: ...

    def __get_repr(self) -> str: ...

    def __get_str(self) -> str: ...

class _Sample:
    def __init__(self, sample: Table, gc: Optional[Cluster] = None) -> None: ...

    def __str__(self) -> str: ...

    def __repr__(self) -> str: ...

    def __get_repr(self) -> str: ...

    def to_pandas(
        self, overwrite: bool = False, *args, **kwargs
    ) -> "pandas.DataFrame": ...

    def to_table(self, *args) -> Table: ...
