from builtins import object

class Potential(object):
  """Class used to describe a potential to the :func:`.writePotentials()` function.

  Potential objects encapsulate a python function or callable which is used by
  the :meth:`.energy` method to calculate potential energy. This callable is also
  used when calculating forces :math:`\\frac{-dU}{dr}` through the  :meth:`.force` method.
  Forces are calculated by using a simple finite difference method to
  find the linear gradient of the energy for a given separation."""

  def __init__(self, speciesA, speciesB, potentialFunction):
    """Create a Potential object from a python function or callable that returns energy at a given separation.

    :param speciesA: Label of first species in the potential pair
    :type spciesA: str
    :param speciesB: Label of second species in the potential pair
    :type speciesB: str
    :param potentialFunction: Python callable which accepts a single parameter (separation) and returns energy for that separation"""
    self._speciesA = speciesA
    self._speciesB = speciesB
    self._potentialFunction = potentialFunction

  @property
  def speciesA(self):
    return self._speciesA

  @property
  def speciesB(self):
    return self._speciesB

  @property
  def potentialFunction(self):
    return self_potentialFunction

  def energy(self, r):
    """:param r: Separation
       :return: Energy for given separation"""
    return self._potentialFunction(r)

  def force(self, r, h = 0.1e-5):
    """Calculate force for this potential at a given separation.

       Force is calculated as -dU/dr by performing numerical differentiation of function.

       :param r: Separation
       :type r: float
       :param h: Distance increment used when calculating energy derivative centred on r
       :type h: float
       :return: -dU/dr at given separation
       :rtype: float"""
    r1 = r - (h / 2.0)
    r2 = r + (h / 2.0)

    dr = r2 -r1
    dU = self.energy(r2) - self.energy(r1)

    dUdr = -dU / dr
    return dUdr
