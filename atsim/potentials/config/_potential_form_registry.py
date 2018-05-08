import inspect
import itertools

from ._common import PotentialFormTuple
from ._common import PotentialFormSignatureTuple
from ._common import Potential_Form_Registry_Exception
from ._common import ConfigParserMissingSectionException

from ._python_potential_function import _Python_Potential_Function
from ._potential_form import Potential_Form
from ._cexprtk_potential_function import _Cexptrk_Potential_Function
from ._potential_form import Potential_Form

class Potential_Form_Registry(object):
  """Factory class that takes [Potential-Form] definitions
  from ConfigParser and turns them into Potential_Form objects"""

  def __init__(self, cfg, register_standard = False):
    """:param cfg: ConfigParser instance.
       :param register_standard: If `True` then functions contained in atsim.potentials.potentialfunctions
          are registered with this object with the `as.` namespace prefix."""

    self._potential_forms = {}

    if register_standard:
      self._potential_forms.update(self._register_standard())
    
    try:
      definitions = cfg.potential_form
      self._potential_forms.update(self._build_potential_forms(definitions))
      self._register_with_each_other()
      self._definitions = definitions
    except ConfigParserMissingSectionException:
      pass

  def _register_standard(self):
    from .. import potentialfunctions
    potential_forms = {}
    for name, pyfunc in inspect.getmembers(potentialfunctions, inspect.isfunction):
      name = "as."+name
      argspec = inspect.getargspec(pyfunc)
      d = PotentialFormTuple(signature = PotentialFormSignatureTuple(name, argspec.args), expression = "")
      func = _Python_Potential_Function(d, pyfunc)
      pf = Potential_Form(func)
      potential_forms[name] = pf
    return potential_forms

  def _build_potential_forms(self, definitions):
    potential_forms = {}
    for d in definitions:
      if d.signature.label in potential_forms:
        raise Potential_Form_Registry_Exception("Two potential forms have the same label: '{0}'".format(d.signature.label))
      func = _Cexptrk_Potential_Function(d)
      pf = Potential_Form(func)
      potential_forms[d.signature.label] = pf
    return potential_forms

  def _register_with_each_other(self):
    # So that each function can rely on other custom functions, add each function to every other
    # function's symbol table
    pairs = list(itertools.permutations(self._potential_forms.values(), 2))
    # print([(a.signature.label, b.signature.label) for (a,b) in pairs])
    for a,b in pairs:
      a._function.register_function(b._function)

  @property
  def registered(self):
    """Returns the labels for the potentials registered here."""
    return sorted(self._potential_forms.keys())

  def __getitem__(self, k):
    return self._potential_forms[k]