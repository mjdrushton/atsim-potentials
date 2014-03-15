import unittest

import sys

from atsim import potentials
from atsim.potentials import _dlpoly_writeTABLE

class TestWriteTable(unittest.TestCase):

  def testForce(self):
    """Check that dlpoly.writeTABLE._calculateForce() wraps potentials.Potential.force to return -r dU/dr values needed by DL_POLY"""
    potfunc = potentials.buck(1388.773, 2.76, 175)
    pot = potentials.Potential("A", "B", potfunc)
    expect = -66990.72564944
    actual = _dlpoly_writeTABLE._calculateForce(pot, 0.5, h=0.001)
    self.assertAlmostEqual(expect, actual)
