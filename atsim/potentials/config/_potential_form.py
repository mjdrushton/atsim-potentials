from ._common import Potential_Form_Exception
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