
def deriv(r, func, h = 0.1e-5):
  """Evaluates the derivative of a unary callable, `func` at a value of `r`.

  If the object `func` has a unary method `deriv(r)`, this will be used to evauluate
  the derivative (allowing analytical derivatives to be used). 

  If `func` does not have a specific `deriv(r)` method then its numerical-derivative of
  will be taken by calling num_deriv()

  :param r: Value at which derivative of `func` should be evaluated.
  :param func: Function whose derivative is to be evaluated.
  :param h: Step size used when performing numerical differentiation.
  :return: Derivative of func at `r`."""
  if hasattr(func, "deriv"):
    return func.deriv(r)
  return num_deriv(r, func, h)

def num_deriv(r, func, h = 0.1e-5):
  """Returns numerical derivative of the callable `func`

  :param r: Value at which derivative of `func` should be evaluated.
  :param func: Function whose gradient is to be evaluated.
  :param h: Step size used when performing numerical differentiation.
  :return: Numerical derivative of func at `r`."""

  r1 = r-(h/2.0)
  r2 = r+(h/2.0)

  dr = r2 -r1
  dU = func(r2) - func(r1)
  dUdr = dU/dr
  return dUdr

class _GradientWrapper(object):
  """Callable class instantiated by gradient()

  Dynamically adds a deriv() method if wrapped function
  has a deriv2() method"""

  def __init__(self, func, h):
    self._wrapped = func
    self._h = h

    if hasattr(self._wrapped, "deriv2"):
      # Add a deriv() method to this object
      def deriv(self, r):
        return self._wrapped.deriv2(r)
      #... bind this to the object
      self.deriv = deriv.__get__(self)

  def __call__(self, r):
    return deriv(r, self._wrapped, self._h)


def gradient(func, h = 0.1e-5):
  """Function wrapper that returns derivative of func.

  If the callable, `func` provides a `.deriv(r)` method this will be used
  to evaluate the derivative of the function, if not the returned function
  will use num_deriv() in gradient evaluation.

  If the callable additionally provides a `.deriv2(r)` method, representing
  its second derivative, the function returned by this routine will have a `deriv()`
  method which will delegate to func.deriv2() when called.

  By providing .deriv() and .deriv2() on the `func` callable analytical descriptions
  of a potential's first and second derivatives may be specified.

  :param func: Function to be wrapped
  :param h: Step size used when performing numerical differentiation
  :return: Function that returns derivative of func"""
  wrapped = _GradientWrapper(func, h)
  return wrapped

import functools
class _rpartial(functools.partial):
    def __call__(self, *args, **kwargs):
        kw = self.keywords.copy()
        kw.update(kwargs)
        return self.func(*(args + self.args), **kwargs)
        