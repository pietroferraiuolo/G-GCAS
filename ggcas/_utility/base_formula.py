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
        self._variables:List[sp.Basic] = []
        self._formula:sp.Basic = None

    @property
    def formula(self) -> sp.Basic:
        """Return the formula"""
        return self._formula

    @property
    def variables(self) -> List[sp.Basic]:
        """Return the variables"""
        return self._variables

    @abstractmethod
    def compute(self, values:List[ArrayLike]) -> ArrayLike:
        """Compute the formula"""
        pass

    @abstractmethod
    def _get_formula(self) -> sp.Basic:
        """Return the formula"""
        pass