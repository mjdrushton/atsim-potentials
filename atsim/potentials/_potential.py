from builtins import object

from ._util import gradient

class Potential(object):
  """Class used to describe a potential to the :func:`.writePotentials()` function.

  Potential objects encapsulate a python function or callable which is used by
  the :meth:`.energy` method to calculate potential energy. This callable is also
  used when calculating forces :math:`\\frac{-dU}{dr}` through the  :meth:`.force` method.
  Forces are calculated by using a simple finite difference method to
  find the linear gradient of the energy for a given separation."""

  def __init__(self, speciesA, speciesB, potentialFunction, h = 1e-6):
    """Create a Potential object from a python function or callable that returns energy at a given separation.

    :param speciesA: Label of first species in the potential pair
    :type spciesA: str
    :param speciesB: Label of second species in the potential pair
    :type speciesB: str
    :param potentialFunction: Python callable which accepts a single parameter (separation) and returns energy for that separation.
    :param h: Distance increment used when calculating numerical derivative of energy to calculate forces in .force() method 
             (if potentialFunction doesn't supply analytical derivative through it's .deriv() method)."""
    self._speciesA = speciesA
    self._speciesB = speciesB
    self._potentialFunction = potentialFunction
    self._derivFunction = gradient(self._potentialFunction, h)

  @property
  def speciesA(self):
    return self._speciesA

  @property
  def speciesB(self):
    return self._speciesB

  @property
  def potentialFunction(self):
    return self._potentialFunction

  def energy(self, r):
    """:param r: Separation
       :return: Energy for given separation"""
    return self._potentialFunction(r)

  def force(self, r):
    """Calculate force for this potential at a given separation.

    If this object's potentialFunction has a .deriv() method this will be used to calculate force (allowing analytical derivatives
    to be specified).

    If potentialFunction doesn't have a deriv method then a numerical derivative of the potential function will be returned instead.

    :param r: Separation
    :type r: float
       
    :return: -dU/dr at given separation
    :rtype: float"""
    return -self._derivFunction(r)
