
class _Python_Potential_Function(object):
  """Callable that can be added to cexprtk symbol_table"""

  def __init__(self, potential_form_tuple, pyfunc):
    self._potential_form_tuple = potential_form_tuple
    self._pyfunc = pyfunc

  def register_function(self, func):
    pass
  
  def __call__(self, *args):
    return self._pyfunc(*args)
