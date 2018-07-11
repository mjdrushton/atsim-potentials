from ._common import Potential_Form_Exception

class Potential_Form(object):
  
  def __init__(self, potential_function):
    """Create Potential_Form object.

    :param potential_function: _Potential_Function object."""
    self.potential_definition = potential_function._potential_form_tuple
    self._function = potential_function

  @property
  def signature(self):
    return self.potential_definition.signature

  @property
  def expression(self):
    return self.potential_definition.expression

  def __call__(self, *args):
    if not len(args) == len(self.signature.parameter_names)-1:
      raise Potential_Form_Exception(
        "Potential form requires {required} arguments but {actual} were provided for potential: {name} {args}".format(
          required = len(self.signature.parameter_names)-1, 
          actual = len(args),
          name = self.signature.label,
          args = " ".join([str(a) for a in args])))

    def f(r):
      funcargs = [r]
      funcargs.extend(args)
      return self._function(*funcargs)
    return f