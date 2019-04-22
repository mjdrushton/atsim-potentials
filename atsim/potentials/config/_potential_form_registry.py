import inspect
import itertools

from ._common import PotentialFormTuple
from ._common import PotentialFormSignatureTuple
from ._common import Potential_Form_Registry_Exception
from ._common import ConfigParserMissingSectionException
from ._common import make_potential_form_tuple_from_function

from ._table_form_builder import Table_Form_Builder
from ._python_potential_function import _Python_Potential_Function
from ._potential_form import Potential_Form, Existing_Potential_Form
from ._cexprtk_potential_function import _Cexptrk_Potential_Function
from ._potential_form import Potential_Form

from ..potentialforms import _iscallable


class Potential_Form_Registry(object):
  """Factory class that takes [Potential-Form] and [Table-Form] definitions
  from ConfigParser and turns them into Potential_Form objects"""

  _standard_namespace = "as."

  def __init__(self, cfg, register_standard = False, register_pymath_functions = False):
    """:param cfg: ConfigParser instance.
       :param register_standard: If `True` then functions contained in atsim.potentials.potentialfunctions
          are registered with this object with the `as.` namespace prefix.
       :param register_pymath_functions: If `True` make functions from the python math module available in cexprtk expressions."""

    self._potential_forms = {}

    if register_standard:
      self._potential_forms.update(self._register_standard())

    self._potential_forms.update(self._build_table_forms(cfg.table_form))

    try:
      definitions = cfg.potential_form
      self._potential_forms.update(self._build_potential_forms(definitions))
      self._definitions = definitions
    except ConfigParserMissingSectionException:
      definitions = []

    if register_pymath_functions:
      self._register_pymath_functions()

    self._register_with_each_other()

    self._definitions = definitions

    if register_standard:
      self._register_from_potentialforms(self._potential_forms)

  def _make_standard_name(self, name):
    return self._standard_namespace + name

  def _register_standard(self):
    from .. import potentialfunctions
    potential_forms = {}
    for name, pyfunc in inspect.getmembers(potentialfunctions, _iscallable):
      name = self._make_standard_name(name)
      d = make_potential_form_tuple_from_function(name, pyfunc)
      func = _Python_Potential_Function(d, pyfunc)
      pf = Potential_Form(func)
      potential_forms[name] = pf
    return potential_forms

  def _register_from_potentialforms(self, potential_forms):
    from .. import potentialforms
    for name, potential_form in inspect.getmembers(potentialforms, _iscallable):
      name = self._make_standard_name(name)
      if not name in potential_forms:
        pf = Existing_Potential_Form(name, potential_form)
        potential_forms[name] = pf

  def _build_potential_forms(self, definitions):
    potential_forms = {}
    for d in definitions:
      if d.signature.label in potential_forms:
        raise Potential_Form_Registry_Exception("Two potential forms have the same label in [Potential-Form] section: '{0}'".format(d.signature.label))
      func = _Cexptrk_Potential_Function(d)
      pf = Potential_Form(func)
      potential_forms[d.signature.label] = pf
    return potential_forms

  def _build_table_forms(self, definitions):
    table_forms = {}

    builder = Table_Form_Builder()

    for d in definitions:
      if d.name in self._potential_forms:
        raise Potential_Form_Registry_Exception("Two potential forms have the same label in [Potential-Form] section: '{0}'".format(d.signature.label))

      pf = builder.create_potential_form(d)
      table_forms[d.name] = pf
    return table_forms


  def _register_with_each_other(self):
    # So that each function can rely on other custom functions, add each function to every other
    # function's symbol table
    pairs = list(itertools.permutations(self._potential_forms.values(), 2))
    for a,b in pairs:
      a.potential_function.register_function(b.potential_function)

  def _register_pymath_functions(self):
    # mathfuncs = ["factorial"]
    from . import _pymath

    new_mathfuncs = []
    # import pdb;pdb.set_trace()
    namespace = "pymath"

    for name, pyfunc in inspect.getmembers(_pymath, inspect.isfunction):
      if name.startswith("_"):
        continue
      label = "{}.{}".format(namespace, name)
      d = make_potential_form_tuple_from_function(label, pyfunc)
      func = _Python_Potential_Function(d, pyfunc)
      new_mathfuncs.append(func)

    for pform in self._potential_forms.values():
      for pyfunc in new_mathfuncs:
        pform.potential_function.register_function(pyfunc)


  @property
  def registered(self):
    """Returns the labels for the potentials registered here."""
    return sorted(self._potential_forms.keys())

  def __getitem__(self, k):
    return self._potential_forms[k]
