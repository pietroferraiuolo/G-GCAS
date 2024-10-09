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

class BaseFormula(ABC):
    """
    Base class for formulas calsses
    """
    def __init__(self):
        """The constructor"""
        self._variables:List[sb]    = []
        self._formula:sb            = None
        self._errFormula:sb         = None
        self._errVariables:List[sb] = None
        self._values:ArrayLike      = None
        self._errors:ArrayLike      = None

    @property
    def formula(self) -> sb:
        """Return the formula"""
        return self._formula

    @property
    def variables(self) -> List[sb]:
        """Return the variables"""
        return self._variables

    @property
    def values(self) -> ArrayLike:
        """Return the values"""
        return "Not computed" if self._values is None else self._values
    
    @property
    def errors(self) -> ArrayLike:
        """Return the errors"""
        return "Not computed" if self._errors is None else self._errors
    
    @property
    def errFormula(self) -> sb:
        """Return the error formula"""
        return self._errFormula

    @property
    def errVariables(self) -> List[sb]:
        """Return the error variables"""
        return self._errVariables

    @abstractmethod
    def compute(self, values:List[ArrayLike]) -> ArrayLike:
        """Compute the formula"""
        pass

    @abstractmethod
    def compute_error(self, values:List[ArrayLike], errors:List[ArrayLike]) -> ArrayLike:
        """Compute the error of the formula"""
        pass

    @abstractmethod
    def _get_formula(self) -> sb:
        """Return the formula"""
        pass