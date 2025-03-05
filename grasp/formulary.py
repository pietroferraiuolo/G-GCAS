import os
import sympy as _sp
from sympy.parsing import latex as _latex, sympy_parser as _symparser
from grasp._utility.base_formulary import BaseFormulary as _BaseFormulary
from grasp._utility.folder_paths import (
    SYS_DATA_FOLDER as _sdf, 
    FORMULARY_BASE_FILE as _fbf
)


# logical check for variables presence in the class
# list(self.formulas['Angular Separation'].free_symbols)[-1].name in [s.name for s in self.formulas['Angular Separation'].free_symbols]


class Formulary(_BaseFormulary):

    def __init__(
        self,
        name: str = "",
        formula_names: list = [],
        formulas: list = []
    ):
        """The constructor"""
        self.name = name
        self.symbols = set({})
        self.functions = set({})
        self.formulas = dict(zip(formula_names, formulas))
        super().__init__()


    def display_all(self):
        """
        Display all formulas in the current formulary instance.
        """
        for name, formula in self.formulas.items():
            if isinstance(formula, _sp.Equality):
                lhs, rhs = formula.args
                print(f"\n{name}\n{lhs} = {rhs}")
            else:
                print(f"\n{name}\n{formula}")


    def load_formulary(self, filename: str = _fbf):
        """
        Load a formulary from a file.
        
        Parameters
        ----------
        filename : str
            The name of the file to load the formulary from.
        
        """
        if not '.frm' in filename:
            filename += '.frm'
        if len(filename.split('/')) == 1:
            filename = os.path.join(_sdf, filename)
        self.name = filename.split('/')[-1].split('.')[0]
        with open(filename, "r") as frm:
            content = frm.readlines()
        self._type = content[0].split(":")[1].strip()
        for i in range(1, len(content), 2):
            if content[i] in ["\n", ''] or '#' in content[i]:
                continue
            name = content[i].strip()
            formula = content[i + 1].strip()
            if self._type == 'latex':
                self.formulas[name] = _latex.parse_latex(formula)
            elif self._type == 'sympy':
                self.formulas[name] = _symparser.parse_expr(formula)
            self.symbols.update(self.formulas[name].free_symbols)
            self.functions.update(self.formulas[name].atoms(_sp.Function))


    def save_formulary(self, filename: str):
        """
        Save the formulary to a file.
        
        Parameters
        ----------
        filename : str
            The name of the file to save the formulary to.
        
        """
        if not '.frm' in filename:
            filename += '.frm'
        if len(filename.split('/')) == 1:
            filename = os.path.join(_sdf, filename)
        with open(filename, "w") as frm:
            frm.write("type: sympy\n")
            for name, formula in self.formulas.items():
                frm.write(f"{name}\n{str(formula)}\n")


def load_formulary(filename: str = _fbf) -> Formulary:
    """"""
    f = Formulary()
    f.load_formulary(filename)
    return f
