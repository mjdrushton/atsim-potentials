import sys

from atsim import potentials
from atsim.potentials import _dlpoly_writeTABLE

import pytest

def testForce():
  """Check that dlpoly.writeTABLE._calculateForce() wraps potentials.Potential.force to return -r dU/dr values needed by DL_POLY"""
  potfunc = potentials.buck(1388.773, 2.76, 175)
  pot = potentials.Potential("A", "B", potfunc)
  expect = pytest.approx(-66990.0984477696)
  actual = _dlpoly_writeTABLE._calculateForce(pot, 0.5)
  assert expect == actual
