from __future__ import division

def gradient(func, h = 0.1e-5):
  """Function wrapper that returns derivative of func.

  @param func Function to be wrapped
  @param h Step size used when performing numerical differentiation
  @return Function that returns derivative of func"""
  def wrapped(r):
    r1 = r-(h/2.0)
    r2 = r+(h/2.0)

    dr = r2 -r1
    dU = func(r2) - func(r1)
    dUdr = dU/dr
    return dUdr
  return wrapped
