from atsim import potentials

import math

import pytest

def test_potclass_testEnergy():
  """Check potentials.Potential.energy"""
  potfunc = potentials.buck(1388.773, 2.76, 175)
  pot = potentials.Potential("A", "B", potfunc)
  assert pytest.approx(-10041.34343169, abs = 1e-5) == pot.energy(0.5)

def test_potclass_testForce_numderiv():
  """Check potentials.Potential.force() method"""
  def energyfunc(r):
    return 1388.773 * math.exp( -r/0.362318841) -  175.0/r**6
  
  pot = potentials.Potential("A", "B", energyfunc, h= 0.1e-5)
  expect = 7.151
  actual = pot.force(2.0)
  assert pytest.approx(expect, abs = 1e-3) == actual

def test_potclass_testForceDeriv():
  """Check that analytical derivatives are used when available"""
  
  # Create a callable that returns an analytical derivative which is very difficult to the (correct)
  # numerical derivative - then look for this when the force() metho of the potential object is called.
  
  deriv_force = 24.0

  class AnalyticalDerivative(object):

    def __call__(self, r):
      return -3.0 * r + 5.0

    def deriv(self, r):
      return deriv_force

  potfunc = AnalyticalDerivative()

  pot = potentials.Potential("A", "B", potfunc)
  expect = 2.0
  actual = pot.energy(1.0)
  assert pytest.approx(expect) == actual

  # Now check the force.
  assert pytest.approx(-deriv_force) == pot.force(1.0)







