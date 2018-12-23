from ._common import Potential_Form_Exception
from ..potentialforms import _FunctionFactory


class Potential_Form(object):
  
  def __init__(self, potential_function):
    """Create Potential_Form object.

    :param potential_function: _Potential_Function object."""
    self.potential_definition = potential_function._potential_form_tuple
    self._potential_function = potential_function
    self._functionfactory = _FunctionFactory(potential_function)

  @property
  def potential_function(self):
    return self._potential_function

  @property
  def signature(self):
    return self.potential_definition.signature

  @property
  def expression(self):
    return self.potential_definition.expression

  def __call__(self, *args):
    if not self.signature.is_varargs and not len(args) == len(self.signature.parameter_names)-1:
      raise Potential_Form_Exception(
        "Potential form requires {required} arguments but {actual} were provided for potential: {name} {args}".format(
          required = len(self.signature.parameter_names)-1, 
          actual = len(args),
          name = self.signature.label,
          args = " ".join([str(a) for a in args])))

    f = self._functionfactory(*args)
    return f