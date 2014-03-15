import unittest
import StringIO
import os

from atsim import potentials

import test_lammpsWriteTABLE
import test_lammpsWriteEAM
import test_dlpoly_writeTABEAM
import test_dlpoly_writeTABLE

class TestPotentialFunctions(unittest.TestCase):
  """Tests for the potential function factories defined in atsim.potentials"""

  def testBuck(self):
    r = 0.5
    a = 1388.773
    rho = 2.76
    c = 175.0

    potfunc = potentials.buck(a,rho,c)
    self.assertAlmostEqual(-10041.34343, potfunc(r), places=5)

  def testHbnd(self):
    a = 300.0
    b = 20.0

    potfunc = potentials.hbnd(a,b)
    self.assertAlmostEqual(1208320, potfunc(0.5))


  def testPlus(self):
    potfunc = potentials.plus( potentials.buck(1388.773, 2.76, 175.0), potentials.hbnd(300.0, 20.0))
    self.assertAlmostEqual(1198278.65656831, potfunc(0.5))

class TestPotentialClass(unittest.TestCase):
  def testEnergy(self):
    """Check potentials.Potential.energy"""
    potfunc = potentials.buck(1388.773, 2.76, 175)
    pot = potentials.Potential("A", "B", potfunc)
    self.assertAlmostEquals(-10041.34343169, pot.energy(0.5), places = 5)

  def testForce(self):
    """Check potentials.Potential.force() method"""
    energyfunc = potentials.buck(1388.773, 0.362318841, 175.0)
    pot = potentials.Potential("A", "B", energyfunc)
    expect = 7.151
    actual = pot.force(2.0, h= 0.1e-5)
    self.assertAlmostEquals(expect, actual, places = 3)

class TestWritePotentials(unittest.TestCase):
  """Tests for the atsim.potentials.writePotentials"""

  def setUp(self):
    f_OO = potentials.potentialforms.buck(1633.00510, 0.327022, 3.948790)
    f_UU = potentials.potentialforms.buck(294.640000, 0.327022, 0.0)
    potential_objects = [
        potentials.Potential('O', 'O', f_OO),
        potentials.Potential('U', 'U', f_UU)]
    self.potential_objects = potential_objects

  def testDL_POLY(self):
    """Test DL_POLY tabulation"""
    sio = StringIO.StringIO()
    potentials.writePotentials(
       'DL_POLY',
       self.potential_objects,
       6.5, 6500,
       out = sio)
    sio.seek(0)
    sio.next()
    line = sio.next()
    tokens = line.split()
    delpot, cutpot, ngrid = tokens

    self.assertAlmostEquals(6.5/float(6500), float(delpot), places = 5)
    self.assertAlmostEquals(6.5, float(cutpot))
    self.assertEquals(6500, int(ngrid))

  def testLAMMPS(self):
    """Test potentials.writePotentials() for LAMMPS"""

    pota = potentials.Potential("A", "B", lambda x: x)
    potb = potentials.Potential("C", "D", lambda x: 6.0-x)

    expect=["A-B",
            "N 6 R 1.00000000 6.00000000",
            "",
            "1 1.00000000 1.00000000 -1.00000000",
            "2 2.00000000 2.00000000 -1.00000000",
            "3 3.00000000 3.00000000 -1.00000000",
            "4 4.00000000 4.00000000 -1.00000000",
            "5 5.00000000 5.00000000 -1.00000000",
            "6 6.00000000 6.00000000 -1.00000000",
            "",
            "C-D",
            "N 6 R 1.00000000 6.00000000",
            "",
            "1 1.00000000 5.00000000 1.00000000",
            "2 2.00000000 4.00000000 1.00000000",
            "3 3.00000000 3.00000000 1.00000000",
            "4 4.00000000 2.00000000 1.00000000",
            "5 5.00000000 1.00000000 1.00000000",
            "6 6.00000000 0.00000000 1.00000000"]

    sbuild = StringIO.StringIO()
    potentials.writePotentials('LAMMPS', [pota,potb], 6.0, 6, sbuild)
    sbuild.seek(0)
    actual = sbuild.readlines()
    msg = "%s != %s" % (expect, actual)

    self.assertEquals(len(expect), len(actual), msg = msg)
    for e,a in zip(expect, actual):
      self.assertEquals(os.linesep, a[-1])
      a = a[:-1]
      self.assertEquals(e,a)
