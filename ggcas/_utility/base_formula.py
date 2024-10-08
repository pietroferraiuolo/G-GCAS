"""
Author(s)
---------
    - Pietro Ferraiuolo : Written in 2024

Description
-----------
Base class for formulas calsses, used in the 'ggcas.functions' module.
"""
from abc import ABC, abstractmethod
from typing import List, Dict, Any, Union
from numpy.typing import ArrayLike
import sympy as sp

class BaseFormula(ABC):
    """
    Base class for formulas calsses
    """
    def __init__(self):
        """The constructor"""
        self._variables:List[sp.Basic]  = []
        self._formula:sp.Basic          = None
        self._errFormula:sp.Basic       = None
        self._errFormula:sp.Basic       = None
        self._values:ArrayLike          = None
        self._errors:ArrayLike          = None

    @property
    def formula(self) -> sp.Basic:
        """Return the formula"""
        return self._formula

    @property
    def variables(self) -> List[sp.Basic]:
        """Return the variables"""
        return self._variables

    @property
    def values(self) -> ArrayLike:
        """Return the values"""
        return self._values
    
    @property
    def errors(self) -> ArrayLike:
        """Return the errors"""
        return self._errors
    
    @property
    def errFormula(self) -> sp.Basic:
        """Return the error formula"""
        return self._errFormula

    @abstractmethod
    def compute(self, values:List[ArrayLike]) -> ArrayLike:
        """Compute the formula"""
        pass

    @abstractmethod
    def compute_error(self, values:List[ArrayLike], errors:List[ArrayLike]) -> ArrayLike:
        """Compute the error of the formula"""
        pass

    @abstractmethod
    def _get_formula(self) -> sp.Basic:
        """Return the formula"""
        pass