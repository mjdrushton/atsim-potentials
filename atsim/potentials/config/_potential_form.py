from ._common import Potential_Form_Exception
from ._common import make_potential_form_tuple_from_function

from ..potentialforms import _FunctionFactory


class _Check_Call(object):
  """Class to check arguments for calls to functions stored in Potential_Form"""

  def __init__(self, signature, is_func_call = False):
    self.signature = signature
    self.is_func_call = is_func_call

  def args_valid(self, *args):
    return self.signature.is_varargs or len(args) == self.required_arg_len()

  def required_arg_len(self):
    argl = len(self.signature.parameter_names)
    if not self.is_func_call:
      argl = argl-1
    return argl

  def recommended_usage(self):
    if self.is_func_call:
      argstring = ",".join(list(self.signature.parameter_names))
      argstring = "({})".format(argstring)
    else:
      argstring = " ".join(list(self.signature.parameter_names)[1:])
      argstring = " "+argstring
    recommendation = "expected usage: '{name}{args}'".format(name = self.signature.label, args = argstring)
    return recommendation

  def how_used(self, *args):
    if self.is_func_call:
      return "("+",".join([str(a) for a in args])+")"
    else:
      return " "+" ".join([str(a) for a in args])

  def __call__(self, *args):
    if not self.args_valid(*args):
      form_or_function = "form"
      if self.is_func_call:
        form_or_function = "function"
      msg = "Potential {form} requires {required} arguments but {actual} were provided for '{name}{args}'".format(
          form = form_or_function,
          required = self.required_arg_len(), 
          actual = len(args),
          name = self.signature.label,
          args = self.how_used(*args))
      recommendation = self.recommended_usage()
      msg = "{} ({})".format(msg, recommendation)
      raise Potential_Form_Exception(msg)


class Potential_Form(object):
  """Wraps cexprtk and python functions so that they can be used
  as potential-form style function factories in Potential_Form_Registry"""
  
  def __init__(self, potential_function):
    """Create Potential_Form object.

    :param potential_function: _Potential_Function object."""
    self.potential_definition = potential_function._potential_form_tuple
    self._potential_function = potential_function
    self._functionfactory = _FunctionFactory(potential_function)
    self._check_call = _Check_Call(self.signature)


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
    self._check_call(*args)
    f = self._functionfactory(*args)
    return f

class Existing_Potential_Form(object):
  """Wraps an object from the atsim.potentials.potentialforms module so that 
  it can be registered directly with Potential_Form_Registry and so used from
  config files. 
  
  This is the case for objects that are registered in potentialforms
  but not in the potentialfunctions module"""


  def __init__(self, name, potential_form):
    """:param name: String giving potential form name (without namespace prexfix)
       :param potential_form: Function factory of the kind provided in potentialforms module"""
    self._potential_form = potential_form
    self.potential_definition = self._make_definition(name)
    self._check_call = _Check_Call(self.signature)

  def _make_definition(self, name):
    pft = make_potential_form_tuple_from_function(name, self._potential_form)
    orig_params = pft.signature.parameter_names
    withr_params = ["r"]
    withr_params.extend(orig_params)
    signature = pft.signature._replace(parameter_names = withr_params)
    pft = pft._replace(signature = signature)
    return pft

  @property
  def signature(self):
    return self.potential_definition.signature

  def __call__(self, *args):
    self._check_call(*args)
    f = self._potential_form(*args)
    return f