from ggcas._cluster import Cluster
from ggcas._utility.sample import Sample
from ggcas._utility.base_formula import BaseFormula
from ggcas._query import GaiaQuery
from ggcas.analyzers._Rcode import r2py_models
from typing import TypeVar, Protocol

_BaseFormulaType = TypeVar('_BaseFormulaType', bound=BaseFormula)
class _SupportsBaseFormula(Protocol):
    def __base_formula__(self) -> BaseFormula: ...

GSample = Sample
GCluster = Cluster
GFormula = _BaseFormulaType
GQuery = GaiaQuery
GRRegression = r2py_models.RegressionModel
GRRModel = r2py_models.GMModel
