"""Functions representing different potential forms.

The functions contained herein are function factories returning
a function that takes separation as its sole argument.

See :ref:`list-of-potential-forms` for descriptions of these potential forms.

"""

from . import potentialfunctions
import inspect
import sys

from ._util import _rpartial

try:
  from collections.abc import Callable
except ImportError:
  from collections import Callable

def potential(func):
  """Decorator for callables that should be tagged as potential-forms or potential-functions"""
  func.is_potential = True
  return func

def is_potential(obj):
  """Identifies if an object is a potential-form or potential-function"""
  return hasattr(obj, "is_potential") and obj.is_potential

def _iscallable(obj):
  return (not obj is Callable) and isinstance(obj, Callable) and (not inspect.isclass(obj)) and is_potential(obj)

class _FunctionFactory(object):

  is_potential = True

  def __init__(self, func):
    self._func = func

  def __call__(self, *args):
    wrapper = _rpartial(self._func, *args)
    if hasattr(self._func, "deriv"):
      wrapper.deriv = _rpartial(self._func.deriv, *args)
    if hasattr(self._func, "deriv2"):
      wrapper.deriv2 = _rpartial(self._func.deriv2, *args)
    return wrapper

@potential
def buck4(A, rho, C, r_detach, r_min, r_attach):
  """Returns a potential form describing the four-range Buckingham potential.

  The potential form is:

  .. math::

    V(r_{ij}) = 
    \\begin{cases}
      A \\exp(-r_{ij}/\\rho)                                                 , & 0 \\leq r_{ij} \\leq  r_\\text{detach}\\\\
      a_0 + a_1 r_{ij} +a_2 r_{ij}^2+a_3 r_{ij}^3+a_4 r_{ij}^4+a_5 r_{ij}^5, & r_\\text{detach} < r_{ij} < r_\\text{min}\\\\
      b_0 +b_1 *r_{ij}+b_2*r_{ij}^2+b_3*r_{ij}^3                           , & r_\\text{min} \\leq r_{ij} < r_\\text{attach}\\\\
      -\\frac{C}{r_{ij}^6}                                                  , & r_{ij} \\geq r_\\text{attach}\\\\
    \end{cases}

  In other words this is a Buckingham potential in which the Born-Mayer component acts at small separations and the disprsion
  term acts at larger separation. These two parts are linked by a fifth then third order polynomial (with a minimum formed in the spline
  at :math:`r_\text{min}`).

  The spline parameters are subject to the constraints that :math:`V(r_{ij})`, first and second derivatives must be equal at the boundary points
  and the function must have a stationary point at `r_min`.

  .. seealso::

    * :class:`atsim.potentials.Buck4_Spline`
    * :class:`atsim.potentials.Buck4_SplinePotential`

  .. note::

    Due to the complexity of calculating the spline-coefficients this potential form does not have an equivalent in the atsim.potentials.potentialfunctions
    module.


  :param A: A potential parameter.
  :param rho: potential parameter.
  :param C: C parameter.
  :param r_detach: Separation where spline starts.
  :param r_min: Location of stationary point.
  :param r_attach: End of splined region.

  :return: Splined potential."""
  bm = bornmayer(A,rho)
  disp = buck(0.0, 1.0, C)

  from .spline import Buck4_SplinePotential
  pot = Buck4_SplinePotential(bm, disp, r_detach, r_attach, r_min)

  return pot

buck = _FunctionFactory(potentialfunctions.buck)
bornmayer = _FunctionFactory(potentialfunctions.bornmayer)
coul = _FunctionFactory(potentialfunctions.coul)
constant = _FunctionFactory(potentialfunctions.constant)
exponential = _FunctionFactory(potentialfunctions.exponential)
hbnd = _FunctionFactory(potentialfunctions.hbnd)
lj = _FunctionFactory(potentialfunctions.lj)
morse = _FunctionFactory(potentialfunctions.morse)
polynomial = _FunctionFactory(potentialfunctions.polynomial)
sqrt = _FunctionFactory(potentialfunctions.sqrt)
tang_toennies = _FunctionFactory(potentialfunctions.tang_toennies)
zbl = _FunctionFactory(potentialfunctions.zbl)
zero = _FunctionFactory(potentialfunctions.zero)
exp_spline = _FunctionFactory(potentialfunctions.exp_spline)