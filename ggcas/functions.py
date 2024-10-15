"""
Module: functions.py

Author(s)
---------
    - Pietro Ferraiuolo : Written in 2024

Description
-----------
This module provides a set of functions to compute various astronomical and
physical quantities such as angular separation, line-of-sight distance, radial
distances in 2D and 3D, total velocity, and effective gravitational potential.

Functions
---------
- angular_separation(ra0=None, dec0=None)
    Computes the angular separation between two points in the sky.

- los_distance()
    Computes the line-of-sight distance based on parallax.

- radial_distance_2d(analytical_w=False, **params)
    Computes the 2D-projected radial distance of a source from the center of a cluster or given RA/DEC coordinates.

- radial_distance_3d(gc_distance=None, analytical_r2d: bool = False, analytical_w: bool = False)
    Computes the 3D radial distance of a source from the center of a cluster or given RA/DEC coordinates.

- total_velocity()
    Computes the total velocity based on the given velocity components.

- effective_potential(shell: bool = False)
    Computes the effective gravitational potential, optionally considering a shell model.

How to Use
----------
1. Import the module:
    ```python
    import functions
    ```

2. Call the desired function with appropriate parameters:
    ```python
    result = functions.angular_separation(ra0=10.0, dec0=20.0)
    ```

Examples
--------
Example usage of `angular_separation` function:
    ```python
    import functions
    from astropy import units as u

    ra0 = 10.0 * u.deg
    dec0 = 20.0 * u.deg
    result = functions.angular_separation(ra0, dec0)
    print(result)
    ```

Example usage of `los_distance` function:
    ```python
    import functions

    result = functions.los_distance()
    print(result)
    ```

Example usage of `radial_distance_2d` function:
    ```python
    import functions

    result = functions.radial_distance_2d(analytical_w=True, ra0=10.0, dec0=20.0)
    print(result)
    ```

Example usage of `radial_distance_3d` function:
    ```python
    import functions
    from astropy import units as u

    gc_distance = 1000 * u.pc
    result = functions.radial_distance_3d(gc_distance=gc_distance, analytical_r2d=True, analytical_w=True)
    print(result)
    ```

Example usage of `total_velocity` function:
    ```python
    import functions

    result = functions.total_velocity()
    print(result)
    ```

Example usage of `effective_potential` function:
    ```python
    import functions

    result = functions.effective_potential(shell=True)
    print(result)
    ```
"""

import numpy as _np
import sympy as _sp
from typing import List as _List, Union as _Union
from numpy.typing import ArrayLike as _ArrayLike
from astropy import units as _u
from ggcas._utility.base_formula import BaseFormula
from ggcas.analyzers.calculus import (
    compute_numerical_function as _compute_numerical,
    error_propagation as _error_propagation,
    compute_error as _compute_error,
    _data_dict_creation,
)

class AngularSeparation(BaseFormula):

    def __init__(self, ra0: _Union[float, _u.Quantity], dec0: _Union[float, _u.Quantity]):
        """The constructor"""
        super().__init__()
        self.ra0 = ra0
        self.dec0 = dec0
        self._get_formula()

    def _get_formula(self):
        """Analytical formula getter for the angular separation"""
        ra1, dec1 = _sp.symbols("alpha_x \delta_x")
        variables = [ra1, dec1]
        d2r = _np.pi / 180
        if isinstance(self.ra0,_u.Quantity):
            self.ra0 = float(self.ra0 /_u.deg)
        if isinstance(self.dec0,_u.Quantity):
            self.dec0 = float(self.dec0 /_u.deg)
        costerm = _np.cos(self.dec0 * d2r)
        w = (
            2
            * _sp.asin(
                _sp.sqrt(
                    _sp.sin((self.dec0 - dec1) * 0.5 * d2r) ** 2
                    + costerm
                    * _sp.cos(dec1 * d2r)
                    * _sp.sin((self.ra0 - ra1) * 0.5 * d2r) ** 2
                )
            )
            / d2r
        )
        self._formula = w
        self._variables = variables
        error = _error_propagation(
            self._formula, self._variables, correlation=True
        )
        self._errFormula = error["error_formula"]
        errvars = error["error_variables"]
        self._errVariables = errvars['errors'] 
        self._correlations = errvars['corrs']
        return self
    
    def _analitical_formula(self):
        ra0, dec0, ra1, dec1 = _sp.symbols("alpha_0 \delta_0 alpha_x \delta_x")
        formula = (
            2
            * _sp.asin(
                _sp.sqrt(
                    _sp.sin((dec0 - dec1)/2) ** 2
                    + _sp.cos(dec0)
                    * _sp.cos(dec1)
                    * _sp.sin((ra0 - ra1)/2) ** 2
                )
            )
        )
        return formula


class LosDistance(BaseFormula):
    """
    Class for the analytical line-of-sight distance.
    """

    def __init__(self):
        """The constructor"""
        super().__init__()
        self._get_formula()

    def _get_formula(self):
        """Analytical formula getter for the line-of-sight distance"""
        parallax = _sp.symbols("omega")
        r = 1 / parallax
        self._formula = r
        self._variables = [parallax]
        error = _error_propagation(
            self._formula, self._variables, correlation=False
        )
        self._errFormula = error["error_formula"]
        errvars = error["error_variables"]
        self._errVariables = errvars['errors']
        return self
    
    def _analitical_formula(self):
        """
        Return the analytical formula for the line-of-sight distance.
        """
        return self._formula


class RadialDistance2D(BaseFormula):
    """
    Class for the analytical 2D-projected radial distance.
    """

    def __init__(self, gc_distance: _Union[float,_u.Quantity]):
        """The constructor"""
        self.gc_distance = gc_distance
        super().__init__()
        self._get_formula()

    def _get_formula(self):
        """Analytical formula getter for the 2D-projected radial distance"""
        variables = []
        if isinstance(self.gc_distance,_u.Quantity):
            rgc = self.gc_distance.to(_u.pc).value
        else:
            rgc = self.gc_distance
        w = _sp.symbols("vartheta_x0")
        variables.append(w)
        r_2d = rgc * _sp.tan(w)
        self._formula = r_2d
        self._variables = variables
        error = _error_propagation(
            self._formula, self._variables, correlation=True
        )
        self._errFormula = error["error_formula"]
        errvars = error["error_variables"]
        self._errVariables = errvars['errors'] 
        self._correlations = errvars['corrs']
        return self
    
    def _analitical_formula(self):
        """
        Return the analytical formula for the 2D-projected radial distance.
        """
        rgc, w = _sp.symbols("r_gc, vartheta_x0")
        formula = rgc * _sp.tan(w)
        return formula


class RadialDistance3D(BaseFormula):
    """
    Class for the analytical 3D radial distance.
    """

    def __init__(
        self, gc_distance:_Union[float,_u.Quantity]=None
    ):
        """The constructor"""
        super().__init__()
        self.gc_distance = gc_distance
        self.extended = False if gc_distance is None else True
        self._get_formula()

    def _get_formula(self):
        """Analytical formula getter for the 3D radial distance"""
        if self.extended:
            if self.gc_distance is None:
                raise ValueError(
                    "The cluster's distance must be provided if extended=True"
                )
            rx = LosDistance()
            r2d = RadialDistance2D(self.gc_distance)
            variables = rx.variables + r2d.variables
            D = rx.formula - self.gc_distance
            r_3d = _sp.sqrt(D**2 + r2d.formula**2)
        else:
            r2d, d = _sp.symbols("r_2d, D")
            variables = [d, r2d]
            r_3d = _sp.sqrt(d**2 + r2d**2)
        self._formula = r_3d
        self._variables = variables
        error = _error_propagation(
            self._formula, self._variables, correlation=True
        )
        self._errFormula = error["error_formula"]
        errvars = error["error_variables"]
        self._errVariables = errvars['errors'] 
        self._correlations = errvars['corrs']
        return self
    
    def _analitical_formula(self):
        """
        Return the analytical formula for the 3D radial distance.
        """
        r2d, d = _sp.symbols("r_2d, D")
        formula = _sp.sqrt(d**2 + r2d**2)
        return formula


class TotalVelocity(BaseFormula):
    """
    Class for the analytical *cartesian* total velocity.

    For cartesian velocities we intend converted ones from proper motion,
    using a simple geometrical conversion. Also, the total velocity, which is
    :math:`v^2 = v_x^2 + v_y^2 + v_z^2`, is computed as :math:`1.5(v_x^2 + v_y^2)`,
    i.e. in the assumption of isotropic distribution of velocities.
    """

    def __init__(self):
        """The constructor"""
        super().__init__()
        self._get_formula()

    def _get_formula(self):
        """Analytical formula getter for the total velocity"""
        vx, vy = _sp.symbols("v_x, v_y")
        V = (3/2)*(vx**2 + vy**2)
        variables = [vx, vy]
        self._formula = V
        self._variables = variables
        error = _error_propagation(
            self._formula, self._variables, correlation=True
        )
        self._errFormula = error["error_formula"]
        errvars = error["error_variables"]
        self._errVariables = errvars['errors'] 
        self._correlations = errvars['corrs']
        return self
    
    def _analitical_formula(self):
        """
        Return the analytical formula for the total velocity.
        """
        return self._formula


class EffectivePotential(BaseFormula):
    r"""
    Class for the analytical effective gravitational potential.

    The effective potential is defined as :math:`\Phi_{\text{eff}} = \ln B - \ln \left(\frac{\Delta N}{\sqrt{x}}\right) - x`,
    where :math:`\ln B = 16\sqrt{2} A (\pi b r)^2 a (m\sigma)^3 \delta r \delta x`.

    The shell model is also considered, where the effective potential is defined as :math:`\ln B - \ln \left(\frac{\Delta N}{\sqrt{x}}\right) - x`.
    """

    def __init__(self, shell: bool = False):
        """The constructor"""
        super().__init__()
        self.shell = shell
        self._constants = None
        self._get_formula()

    def _get_formula(self):
        """Analytical formula getter for the effective gravitational potential"""
        if self.shell:
            dN, x, lnB = _sp.symbols("\Delta\ N, x, lnB")
            A, a, b, r, m, sigma, dr, dx = _sp.symbols(
                "A, alpha, beta, r_*, m, sigma, dr, dx"
            )
            lnb = (
                16
                * _sp.sqrt(2)
                * A
                * (_sp.pi * b * r) ** 2
                * a
                * (m * sigma) ** 3
                * dr
                * dx
            )
            variables = [dN, x, lnB]
            self._constants = [lnb, A, a, b, r, m, sigma, dr, dx]
            poteff = lnB - _sp.ln(dN / _sp.sqrt(x)) - x
        else:
            x, w = _sp.symbols("x, w")
            variables = [x, w]
            poteff = -_sp.ln(1 - _sp.exp(x - w))
        self._formula = poteff
        self._variables = variables
        self._bkcVariables = self._variables
        error = _error_propagation(
            self._formula, self._variables, correlation=False
        )
        self._errFormula = error["error_formula"]
        errvars = error["error_variables"]
        self._errVariables = errvars['errors']
        return self

    @property
    def constants(self):
        """
        Return the constants used in the effective gravitational potential.
        """
        return self._constants
    
    def set_constants(self, values:_ArrayLike):
        """
        Set the constants used in the effective gravitational potential.

        Parameters
        ----------
        values :_List[ArrayLike]
            The values of the constants to set.
        """
        if self.shell:
            self._variables = self.bkcVariables 
            print(f"WARNING! Be sure that the input data follow this specific order: {self.constants[1:]}")
            c = []
            for x in values:
                c.append([x])
            c = _np.array(c)
            const_dict = _data_dict_creation(self._constants[1:], c)
            lnb = self._constants[0].subs(const_dict)
            self._formula = self._formula.subs(self._variables[2], lnb)
            self._variables.drop(2)
            self._constants = const_dict
            self._constants[_sp.symbols('lnB')] = lnb            
        else:
            print("No constants to set, as shell={self.shell}")
        return self
    
    def _analitical_formula(self):
        """
        Return the analytical formula for the effective gravitational potential.
        """
        return self._formula
    
class CartesianConversion(BaseFormula):
    r"""
    Class for the analytical cartesian conversion.

    The cartesian conversion is defined as :math:`x = \sin(\alpha - \alpha_0) \cos(\delta_0)`,
    :math:`y = \sin(\delta)\cos(\delta_0) - \cos(\delta)\sin(\delta_0)\cos(\alpha - \alpha_0)`,
    and :math:`r = \sqrt{x^2 + y^2}`.
    """

    def __init__(self, ra0=0, dec0=0):
        """The constructor"""
        super().__init__()
        self.ra0 = ra0
        self.dec0 = dec0
        self._get_formula()

    @property
    def x(self):
        """
        Return the x component of the cartesian conversion.
        """
        return self._values[0]
    
    @property
    def y(self):
        """
        Return the y component of the cartesian conversion.
        """
        return self._values[1]
    
    @property
    def r(self):
        """
        Return the r component of the cartesian conversion.
        """
        return self._values[2]

    def _get_formula(self):
        """Analytical formula getter for the cartesian conversion"""
        ra, dec = _sp.symbols("alpha \delta")
        variables = [ra, dec]
        x = _sp.sin(ra - self.ra0) * _sp.cos(self.dec0)
        y = _sp.sin(dec)*_sp.cos(self.dec0) - _sp.cos(dec)*_sp.sin(self.dec0)*_sp.cos(ra - self.ra0)
        r = _sp.sqrt(x**2 + y**2)
        self._formula = [x, y, r]
        self._variables = variables
        return self
    
    def _analitical_formula(self):
        """
        Return the analytical formula for the cartesian conversion.
        """
        ra, dec = _sp.symbols("alpha \delta")
        x = _sp.sin(ra - ra0) * _sp.cos(dec0)
        ra0, dec0 = _sp.symbols("alpha_0 \delta_0")
        y = _sp.sin(dec)*_sp.cos(dec0) - _sp.cos(dec)*_sp.sin(dec0)*_sp.cos(ra - ra0)
        r = _sp.sqrt(x**2 + y**2)
        return [x, y, r]

def cartesian_conversion(data:_List[_ArrayLike], ra0=0, dec0=0) -> _List[_ArrayLike]:
    """
    Convert the proper motion to cartesian velocities.

    Parameters
    ----------
    data :_List[ArrayLike]
        The data to use for the computation.
        The data to provide are:
        - :math:`\alpha*`: the right ascension coordinate.
        - :math:`\delta`: the declination coordinate.
    ra0 : float
        The right ascension of the source.
    dec0 : float
        The declination of the source.

    Returns
    -------
    result :_List[ArrayLike]
        The converted cartesian velocities.
    """
    ra, dec = data
    x = _np.sin(ra - ra0) * _np.cos(dec0)
    y = _np.sin(dec)*_np.cos(dec0) - _np.cos(dec)*_np.sin(dec0)*_np.cos(ra - ra0)
    r = np.sqrt(x**2 + y**2)
    return [x, y, r]


