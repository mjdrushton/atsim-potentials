from ._potential_form import _Check_Call


class _Python_Potential_Function(object):
  """Callable that can be added to cexprtk symbol_table"""

  def __init__(self, potential_form_tuple, pyfunc):
    self._potential_form_tuple = potential_form_tuple
    self._pyfunc = pyfunc
    self._check_call = _Check_Call(potential_form_tuple.signature, True)

    if hasattr(pyfunc, 'deriv'):
      def deriv(self, *args):
        self._check_call(*args)
        return self._pyfunc.deriv(*args)
      self.deriv = deriv.__get__(self)

    if hasattr(pyfunc, 'deriv2'):
      def deriv2(self, *args):
        self._check_call(*args)
        return self._pyfunc.deriv2(*args)
      self.deriv2 = deriv2.__get__(self)

  def register_function(self, func):
    pass
  
  def __call__(self, *args):
    self._check_call(*args)
    return self._pyfunc(*args)

