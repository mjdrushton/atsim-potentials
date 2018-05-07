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
    # Curry all but the first argument
    if not len(args) == len(self.signature.parameter_names)-1:
      raise Potential_Form_Exception(
        "Potential form requires {0} arguments but {1} were provided".format(len(self.signature.parameter_names)-1, len(args)))

    def f(r):
      funcargs = [r]
      funcargs.extend(args)
      return self._function(*funcargs)
    return f