from typing import Optional, Dict, Union, Any
import numpy as np
import sympy as sp
from astropy import units as u

def angular_separation(
    ra0: Optional[Union[float, u.Quantity]] = None,
    dec0: Optional[Union[float, u.Quantity]] = None,
) -> Dict[str, Any]: ...

def los_distance() -> Dict[str, Any]: ...

def radial_distance_2d(analytical_w: bool = False, **params: Any) -> Dict[str, Any]: ...

def radial_distance_3d(
    gc_distance: Optional[Union[float, u.Quantity]] = None,
    analytical_r2d: bool = False,
    analytical_w: bool = False,
) -> Dict[str, Any]: ...

def total_velocity() -> Dict[str, Any]: ...

def effective_potential(shell: bool = False) -> Dict[str, Any]: ...

