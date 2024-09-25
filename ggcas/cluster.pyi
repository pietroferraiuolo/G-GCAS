from typing import Optional, Dict, Any
import astropy.units as u
from astropy.table import Table
import pandas as pd

class Cluster:
    def __init__(self, name: Optional[str] = None, **params: Any) -> None: ...

    def __str__(self) -> str: ...

    def __repr__(self) -> str: ...

    def show_model(self, **kwargs: Any) -> None: ...

    def _load_cluster_parameters(self, name: str) -> pd.Series: ...

    def _load_king_model(self) -> Table: ...

    def __get_repr(self) -> str: ...
    
    def __get_str(self) -> str: ...
