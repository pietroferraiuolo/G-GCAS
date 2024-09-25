from typing import Optional, Union, Tuple, Dict
import numpy as np
import pandas as pd
import ggcas

def scatter_2hist(
    x: np.ndarray, y: np.ndarray, kde: bool = False, **kwargs
) -> None: ...

def colorMagnitude(
    g: Union[float, np.ndarray],
    b_r: Union[float, np.ndarray],
    teff_gspphot: Union[float, np.ndarray],
    **kwargs
) -> None: ...

def properMotion(sample:Union[ggcas._query._Sample, pd.DataFrame], **kwargs) -> None: ...

def spatial(sample:Union[ggcas._query._Sample, pd.DataFrame], **kwargs) -> None: ...

def histogram(
    data: np.ndarray, kde: bool = False, **kwargs
) -> Dict[str, Union[Tuple[np.ndarray, np.ndarray], Tuple[float, float]]]: ...

def scat_xhist(
    x: np.ndarray,
    y: np.ndarray,
    xerr: Optional[Union[float, np.ndarray]] = None,
    **kwargs
) -> Tuple[float, float]: ...

def errorbar(
    data: np.ndarray,
    dataerr: np.ndarray,
    x: Optional[np.ndarray] = None,
    xerr: Optional[np.ndarray] = None,
    **kwargs
) -> None: ...
