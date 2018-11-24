from __future__ import division

def num_gradient(r, func, h = 0.1e-5):
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

def gradient(func, h = 0.1e-5):
  """Function wrapper that returns derivative of func.

  :param func: Function to be wrapped
  :param h: Step size used when performing numerical differentiation
  :return: Function that returns derivative of func"""
  def wrapped(r):
    return num_gradient(r, func, h)
  return wrapped

import functools
class _rpartial(functools.partial):
    def __call__(self, *args, **kwargs):
        kw = self.keywords.copy()
        kw.update(kwargs)
        return self.func(*(args + self.args), **kwargs)
