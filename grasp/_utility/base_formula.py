from __future__ import annotations
"""
Author(s)
---------
    - Pietro Ferraiuolo : Written in 2024

Description
-----------
Base class for formulas calsses, used in the 'grasp.functions' module.
"""

import numpy as np
from typing import List
from sympy import Basic as sb
from numpy.typing import ArrayLike
from abc import ABC, abstractmethod
from grasp.analyzers.calculus import (
    compute_numerical_function as _compute_numerical,
    compute_error as _compute_error,
)


class BaseFormula(ABC):
    """
    Base class for the various formula calsses
    """

    def __init__(self):
        """The constructor"""
        self._variables: List[sb] = []
        self._formula: sb = None
        self._errFormula: sb = None
        self._errVariables: List[sb] = None
        self._correlations: List[sb] = None
        self._values: ArrayLike = None
        self._errors: ArrayLike = None

    def __str__(self) -> str:
        """Return the actual formula as a string"""
        return self._get_str()

    def __repr__(self) -> str:
        """Return the analytical formula as a string"""
        return self._get_str()

    def _get_str(self):
        formula = self._analitical_formula()
        computed = (True if self._values is not None else False) or (
            True if self._errors is not None else False
        )
        return f"{self.__class__.__name__}:\n{formula}\nComputed: {computed}"

    @property
    def formula(self) -> sb:
        """Return the formula"""
        return self._formula

    @property
    def variables(self) -> List[sb]:
        """Return the variables"""
        return self._variables

    @property
    def computed_values(self) -> ArrayLike:
        """Return the values"""
        return "Not computed" if self._values is None else self._values

    @property
    def computed_errors(self) -> ArrayLike:
        """Return the errors"""
        return "Not computed" if self._errors is None else self._errors

    @property
    def error_formula(self) -> sb:
        """Return the error formula"""
        return self._errFormula

    @property
    def error_variables(self) -> List[sb]:
        """Return the error variables"""
        return self._errVariables

    @property
    def correlations(self) -> List[sb]:
        """Return the correlations"""
        return self._correlations

    def compute(
        self,
        data: List[ArrayLike],
        errors: List[ArrayLike] = None,
        correlations: List[ArrayLike] = None,
    ) -> ArrayLike:
        """
        Compute the values of the formula, along with the propagated errors if
        variables errors (and correlation) are provided.

        Parameters
        ----------
        data :_List[ArrayLike]
            The data to use for the computation.
        errors:_List[ArrayLike], optional
            If provided, the propagated errors will be computed.
        correlations:_List[ArrayLike], optional
            If provided, the correlations will be used in the error computation.

        Returns
        -------
        result : ArrayLike
            The computed values.
        """
        print(
            f"""WARNING! Be sure that the input data follow this specific order: 
Data:         {self.variables}"""
        )
        if errors is not None:
            print(
                f"""Errors:       {self._errVariables}
Correlations: {self._correlations}
"""
                + "-" * 30
            )
            variables = self._variables + self._errVariables
            print("Errors:")
            if self._correlations is None:
                self._errors = _compute_error(self._errFormula, variables, data, errors)
            else:
                variables += self._correlations
                if correlations is None:
                    correlations = [
                        np.zeros(len(data[0])) for _ in range(len(self._correlations))
                    ]
                self._errors = _compute_error(
                    self._errFormula, variables, data, errors, corr_values=correlations
                )
        print("Data:")
        self._values = _compute_numerical(self._formula, self._variables, data)
        return self

    @abstractmethod
    def _get_formula(self) -> sb:
        """Return the derived formula of the class"""
        pass

    @abstractmethod
    def _analitical_formula(self) -> str:
        """Return the analitical formula for the class"""
        pass
