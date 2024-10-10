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

import numpy as np
import sympy as sp
from typing import List, Union
from numpy.typing import ArrayLike
from astropy import units as u
from ggcas._utility.base_formula import BaseFormula
from ggcas.analyzers.calculus import (
    compute_numerical_function as _compute_numerical,
    compute_error as _compute_error,
    error_propagation,
)


class AngularSeparation(BaseFormula):

    def __init__(self, ra0: Union[float, u.Quantity], dec0: Union[float, u.Quantity]):
        """The constructor"""
        super().__init__()
        self.ra0 = ra0 if ra0 is not None else sp.symbols("alpha_0")
        self.dec0 = dec0 if dec0 is not None else sp.symbols("\delta_0")
        self._get_formula()

    def _get_formula(self):
        """Analytical formula getter for the angular separation"""
        ra1, dec1 = sp.symbols("alpha_x \delta_x")
        variables = [ra1, dec1]
        d2r = np.pi / 180
        if isinstance(self.ra0, u.Quantity):
            self.ra0 = float(self.ra0 / u.deg)
        if isinstance(self.dec0, u.Quantity):
            self.dec0 = float(self.dec0 / u.deg)
        costerm = np.cos(self.dec0 * d2r)
        w = (
            2
            * sp.asin(
                sp.sqrt(
                    sp.sin((self.dec0 - dec1) * 0.5 * d2r) ** 2
                    + costerm
                    * sp.cos(dec1 * d2r)
                    * sp.sin((self.ra0 - ra1) * 0.5 * d2r) ** 2
                )
            )
            / d2r
        )
        self._formula = w
        self._variables = variables
        error = error_propagation(
            self._formula, self._variables, correlation=True
        )
        self._errFormula = error["error_formula"]
        errvars = error["error_variables"]
        self._errVariables = errvars['errors'] + errvars['correlations']
        return self

    def compute(
            self, 
            data: List[ArrayLike],
            errors:List[ArrayLike]=None,
            correlations:List[ArrayLike]=None
        ) -> ArrayLike:
        f"""
        Compute the angular separation between two points in the sky, along with
        the propagated errors if they are provided.

        Parameters
        ----------
        data : List[ArrayLike]
            The data to use for the computation.
        errors: List[ArrayLike], optional
            If provided, the propagated errors for the {self.__class__.__name__} will be computed.
            Note: the errors can also be computed separately using the `compute_error` method.

        Returns
        -------
        result : ArrayLike
            The computed angular separation.
        """
        if errors is None:
            print(
                f"WARNING! Be sure that the input data follow this specific order: {self.variables}"
            )
            self._values = _compute_numerical(self._formula, self._variables, data)
        else:
            print(
                f"WARNING! Be sure that the input data follow this specific order: {self.errVariables}"
            )
            self._values = _compute_numerical(self._formula, self._variables, data)
            self._errors = _compute_error(
                self._formula, self._variables, data, errors, corr_values=correlations,
            )
        return self

    def compute_error(
        self, data: List[ArrayLike], errors: List[ArrayLike], corr: ArrayLike = None
    ) -> ArrayLike:
        """
        Compute the error of the angular separation.

        Parameters
        ----------
        data : List[ArrayLike]
            The data to use for the computation.
        errors : List[ArrayLike]
            The errors to use for the computation.

        Returns
        -------
        result : ArrayLike
            The computed error of the angular separation.
        """
        print(
            f"WARNING! Be sure that the input data follow this specific order: {self.errVariables}"
        )
        self._errors = _compute_error(
            self._formula, self._variables, data, errors, corr_values=corr
        )
        return self
    
    def _analitical_formula(self):
        ra0, dec0, ra1, dec1 = sp.symbols("alpha_0 \delta_0 alpha_x \delta_x")
        formula = (
            2
            * sp.asin(
                sp.sqrt(
                    sp.sin((dec0 - dec1)/2) ** 2
                    + sp.cos(dec0)
                    * sp.cos(dec1)
                    * sp.sin((ra0 - ra1)/2) ** 2
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
        parallax = sp.symbols("omega")
        r = 1 / parallax
        self._formula = r
        self._variables = [parallax]
        error = error_propagation(
            self._formula, self._variables, correlation=True
        )
        self._errFormula = error["error_formula"]
        errvars = error["error_variables"]
        self._errVariables = errvars['variables'] + errvars['errors'] + errvars['correlations']
        return self

    def compute(self, data: List[ArrayLike]) -> ArrayLike:
        f"""
        Compute the line-of-sight distance based on parallax.

        Parameters
        ----------
        data : List[ArrayLike]
            The data to use for the computation. In this case, the parallax data needs
            to be provided.
        errors: List[ArrayLike], optional
            If provided, the propagated errors for the {self.__class__.__name__} will be computed.
            Note: the errors can also be computed separately using the `compute_error` method.

        Returns
        -------
        result : ArrayLike
            The computed line-of-sight distance.
        """
        print(
            f"WARNING! Be sure that the input data follow this specific order: {self.variables}"
        )
        self._values =  _compute_numerical(self._formula, self._variables, data)
        return self

    def compute_error(
        self, data: List[ArrayLike], errors: List[ArrayLike], corr: ArrayLike = None
    ) -> ArrayLike:
        """
        Compute the error of the line-of-sight distance.

        Parameters
        ----------
        data : List[ArrayLike]
            The data to use for the computation. In this case, the parallax data needs
            to be provided.
        errors : List[ArrayLike]
            The errors to use for the computation.
        corr: ArrayLike, optional
            The correlation values between variables to use for the computation.
            If not provided, the correlation will be set to False.

        Returns
        -------
        result : ArrayLike
            The computed error of the line-of-sight distance.
        """
        print(
            f"WARNING! Be sure that the input data follow this specific order: {self.errVariables}"
        )
        self._errors = _compute_error(
            self._formula, self._variables, data, errors, corr_values=corr
        )
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

    def __init__(self, gc_distance: Union[float, u.Quantity]):
        """The constructor"""
        self.gc_distance = gc_distance
        super().__init__()
        self._get_formula()

    def _get_formula(self):
        """Analytical formula getter for the 2D-projected radial distance"""
        variables = []
        if isinstance(self.gc_distance, u.Quantity):
            rgc = self.gc_distance.to(u.pc).value
        else:
            rgc = self.gc_distance
        w = sp.symbols("vartheta_x0")
        variables.append(w)
        r_2d = rgc * sp.tan(w)
        self._formula = r_2d
        self._variables = variables
        error = error_propagation(
            self._formula, self._variables, correlation=True
        )
        self._errFormula = error["error_formula"]
        errvars = error["error_variables"]
        self._errVariables = errvars['variables'] + errvars['errors'] + errvars['correlations']
        return self

    def compute(self, data: List[ArrayLike]) -> ArrayLike:
        f"""
        Compute the 2D-projected radial distance of a source from the center of a cluster.

        Parameters
        ----------
        data : List[ArrayLike]
            The data to use for the computation.
            The data to provide are:
            - :math:`\theta_x0`: the 2D-projected angular separation from the center of the cluster.
                Can be computed through the 'AngularSeparation' Class.
        errors: List[ArrayLike], optional
            If provided, the propagated errors for the {self.__class__.__name__} will be computed.
            Note: the errors can also be computed separately using the `compute_error` method.

        Returns
        -------
        result : ArrayLike
            The computed 2D-projected radial distance.
        """
        print(
            f"WARNING! Be sure that the input data follow this specific order: {self.variables}"
        )
        self._values = _compute_numerical(self._formula, self._variables, data)
        return self

    def compute_error(
        self, data: List[ArrayLike], errors: List[ArrayLike], corr: ArrayLike = None
    ) -> ArrayLike:
        """
        Compute the error of the 2D-projected radial distance.

        Parameters
        ----------
        data : List[ArrayLike]
            The data to use for the computation.
        errors : List[ArrayLike]
            The errors to use for the computation.

        Returns
        -------
        result : ArrayLike
            The computed error of the 2D-projected radial distance.
        """
        print(
            f"WARNING! Be sure that the input data follow this specific order: {self.errVariables}"
        )
        self._errors = _compute_error(
            self._formula, self._variables, data, errors, corr_values=corr
        )
        return self
    
    def _analitical_formula(self):
        """
        Return the analytical formula for the 2D-projected radial distance.
        """
        rgc, w = sp.symbols("r_gc, vartheta_x0")
        formula = rgc * sp.tan(w)
        return formula


class RadialDistance3D(BaseFormula):
    """
    Class for the analytical 3D radial distance.
    """

    def __init__(
        self, extended: bool = False, gc_distance: Union[float, u.Quantity] = None
    ):
        """The constructor"""
        super().__init__()
        self.gc_distance = gc_distance
        self.extended = extended
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
            r_3d = sp.sqrt(D**2 + r2d.formula**2)
        else:
            r2d, d = sp.symbols("r_2d, D")
            variables = [d, r2d]
            r_3d = sp.sqrt(d**2 + r2d**2)
        self._formula = r_3d
        self._variables = variables
        error = error_propagation(
            self._formula, self._variables, correlation=True
        )
        self._errFormula = error["error_formula"]
        errvars = error["error_variables"]
        self._errVariables = errvars['variables'] + errvars['errors'] + errvars['correlations']
        return self

    def compute(self, data: List[ArrayLike]) -> ArrayLike:
        f"""
        Compute the 3D radial distance of a source from the center of a cluster.

        Parameters
        ----------
        data : List[ArrayLike]
            The data to use for the computation.
        errors: List[ArrayLike], optional
            If provided, the propagated errors for the {self.__class__.__name__} will be computed.
            Note: the errors can also be computed separately using the `compute_error` method.
            
        Returns
        -------
        result : ArrayLike
            The computed 3D radial distance.
        """
        print(
            f"WARNING! Be sure that the input data follow this specific order: {self.variables}"
        )
        self._values = _compute_numerical(self._formula, self._variables, data)
        return self

    def compute_error(
        self, data: List[ArrayLike], errors: List[ArrayLike], corr: ArrayLike = None
    ) -> ArrayLike:
        """
        Compute the error of the 3D radial distance.

        Parameters
        ----------
        data : List[ArrayLike]
            The data to use for the computation.
        errors : List[ArrayLike]
            The errors to use for the computation.

        Returns
        -------
        result : ArrayLike
            The computed error of the 3D radial distance.
        """
        print(
            f"WARNING! Be sure that the input data follow this specific order: {self.errVariables}"
        )
        self._errors = _compute_error(
            self._formula, self._variables, data, errors, corr_values=corr
        )
        return self
    
    def _analitical_formula(self):
        """
        Return the analytical formula for the 3D radial distance.
        """
        r2d, d = sp.symbols("r_2d, D")
        formula = sp.sqrt(d**2 + r2d**2)
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
        vx, vy = sp.symbols("v_x, v_y")
        V = (3/2)*(vx**2 + vy**2)
        variables = [vx, vy]
        self._formula = V
        self._variables = variables
        error = error_propagation(
            self._formula, self._variables, correlation=True
        )
        self._errFormula = error["error_formula"]
        errvars = error["error_variables"]
        self._errVariables = errvars['variables'] + errvars['errors'] + errvars['correlations']
        return self

    def compute(self, data: List[ArrayLike]) -> ArrayLike:
        f"""
        Compute the total velocity based on the given velocity components.

        Parameters
        ----------
        data : List[ArrayLike]
            The data to use for the computation.
        errors: List[ArrayLike], optional
            If provided, the propagated errors for the {self.__class__.__name__} will be computed.
            Note: the errors can also be computed separately using the `compute_error` method.

        Returns
        -------
        result : ArrayLike
            The computed total velocity.
        """
        print(
            f"WARNING! Be sure that the input data follow this specific order: {self.variables}"
        )
        self._values = _compute_numerical(self._formula, self._variables, data)
        return self
    
    def compute_error(
        self, data: List[ArrayLike], errors: List[ArrayLike], corr: ArrayLike = None
    ) -> ArrayLike:
        """
        Compute the error of the total velocity.

        Parameters
        ----------
        data : List[ArrayLike]
            The data to use for the computation.
        errors : List[ArrayLike]
            The errors to use for the computation.

        Returns
        -------
        result : ArrayLike
            The computed error of the total velocity.
        """
        print(
            f"WARNING! Be sure that the input data follow this specific order: {self.errVariables}"
        )
        self._errors = _compute_error(
            self._formula, self._variables, data, errors, corr_values=corr
        )
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
        self._get_formula()

    def _get_formula(self):
        """Analytical formula getter for the effective gravitational potential"""
        if self.shell:
            dN, x, lnB = sp.symbols("\Delta\ N, x, lnB")
            A, a, b, r, m, sigma, dr, dx = sp.symbols(
                "A, alpha, beta, r_*, m, sigma, dr, dx"
            )
            lnb = (
                16
                * sp.sqrt(2)
                * A
                * (sp.pi * b * r) ** 2
                * a
                * (m * sigma) ** 3
                * dr
                * dx
            )
            variables = {
                "x": [dN, x, lnB],
                "const": [lnb, A, a, b, r, m, sigma, dr, dx],
            }
            poteff = lnB - sp.ln(dN / sp.sqrt(x)) - x
        else:
            x, w = sp.symbols("x, w")
            variables = [x, w]
            poteff = -sp.ln(1 - sp.exp(x - w))
        self._formula = poteff
        self._variables = variables
        error = error_propagation(
            self._formula, self._variables, correlation=True
        )
        self._errFormula = error["error_formula"]
        errvars = error["error_variables"]
        self._errVariables = errvars['variables'] + errvars['errors'] + errvars['correlations']
        return self

    def compute(self, data: List[ArrayLike]) -> ArrayLike:
        f"""
        Compute the effective gravitational potential.

        Parameters
        ----------
        data : List[ArrayLike]
            The data to use for the computation.
        errors: List[ArrayLike], optional
            If provided, the propagated errors for the {self.__class__.__name__} will be computed.
            Note: the errors can also be computed separately using the `compute_error` method.

        Returns
        -------
        result : ArrayLike
            The computed effective gravitational potential.

        """
        print(
            f"WARNING! Be sure that the input data follow this specific order: {self.variables}"
        )
        self._values = _compute_numerical(self._formula, self._variables, data)
        return self

    def compute_error(
        self, data: List[ArrayLike], errors: List[ArrayLike], corr: ArrayLike = None
    ) -> ArrayLike:
        """
        Compute the error of the effective gravitational potential.

        Parameters
        ----------
        data : List[ArrayLike]
            The data to use for the computation.
        errors : List[ArrayLike]
            The errors to use for the computation.

        Returns
        -------
        result : ArrayLike
            The computed error of the effective gravitational potential.
        """
        print(
            f"WARNING! Be sure that the input data follow this specific order: {self.errVariables}"
        )
        self._errors = _compute_error(
            self._formula, self._variables, data, errors, corr_values=corr
        )
        return self
    
    def _analitical_formula(self):
        """
        Return the analytical formula for the effective gravitational potential.
        """
        return self._formula
