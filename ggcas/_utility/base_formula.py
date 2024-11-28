"""
Author(s)
---------
    - Pietro Ferraiuolo : Written in 2024

Description
-----------
Base class for formulas calsses, used in the 'ggcas.functions' module.
"""
from abc import ABC, abstractmethod
from typing import List
from numpy.typing import ArrayLike
from sympy import Basic as sb
from ggcas.analyzers.calculus import (
    compute_numerical_function as _compute_numerical,
    compute_error as _compute_error
)

class BaseFormula(ABC):
    """
    Base class for the vrious formula calsses
    """
    def __init__(self):
        """The constructor"""
        self._variables:List[sb]    = []
        self._formula:sb            = None
        self._errFormula:sb         = None
        self._errVariables:List[sb] = None
        self._correlations:List[sb] = None
        self._values:ArrayLike      = None
        self._errors:ArrayLike      = None

    def __str__(self) -> str:
        """Return the actual formula as a string"""
        return self._get_str()
    
    def __repr__(self) -> str:
        """Return the analytical formula as a string"""
        return self._get_str()

    def _get_str(self):
        formula = self._analitical_formula()
        computed = (
            True if self._values is not None else False
        ) or (
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
    
    def compute(self, data:List[ArrayLike], errors:List[ArrayLike]=None, correlations:List[ArrayLike]=None) -> ArrayLike:
        f"""
        Compute the angular separation between two points in the sky, along with
        the propagated errors if they are provided.

        Parameters
        ----------
        data :_List[ArrayLike]
            The data to use for the computation.
        errors:_List[ArrayLike], optional
            If provided, the propagated errors for the {self.__class__.__name__} will be computed.
            Note: the errors can also be computed separately using the `compute_error` method.

        Returns
        -------
        result : ArrayLike
            The computed angular separation.
        """
        if errors is None:
            w_msg = f"""
WARNING! Be sure that the input data follow this specific order: 
Data: {self.variables}
"""
            print(w_msg)
        elif errors is not None:
            if self._correlations is None:
                self.compute_error(data, errors)
            else:
                self.compute_error(data, errors, corr=correlations)
        self._values = _compute_numerical(self._formula, self._variables, data)
        return self

    def compute_error(
        self, data:List[ArrayLike], errors:List[ArrayLike], corr: List[ArrayLike] = None
    ) -> ArrayLike:
        """
        Compute the error of the angular separation.

        Parameters
        ----------
        data :_List[ArrayLike]
            The data to use for the computation.
        errors :_List[ArrayLike]
            The errors to use for the computation.

        Returns
        -------
        result : ArrayLike
            The computed error of the angular separation.
        """
        w_msg = f"""
WARNING! Be sure that the input data follow this specific order: 
Data: {self._variables}
Errors: {self._errVariables}
Correlations: {self._correlations}
"""
        print(w_msg)
        variables = self._variables + self._errVariables
        if corr is not None and self._correlations is not None:
            variables += self._correlations
        self._errors = _compute_error(
            self._errFormula, variables, data, errors, corr_values=corr
        )
        return self

    @abstractmethod
    def _get_formula(self) -> sb:
        """Return the derived formula of the class"""
        pass

    @abstractmethod
    def _analitical_formula(self) -> str:
        """Return the analitical formula for the class"""
        pass