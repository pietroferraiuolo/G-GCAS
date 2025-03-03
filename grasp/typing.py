from grasp._cluster import Cluster
from grasp._utility.sample import Sample
from grasp._utility.base_formula import BaseFormula
from grasp.gaia.query import GaiaQuery
from grasp.analyzers._Rcode import r2py_models
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
