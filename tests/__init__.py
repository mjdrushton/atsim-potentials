import unittest
import StringIO
import os

import atsim_potentials

import test_lammpsWriteTABLE
import test_lammpsWriteEAM
import test_dlpoly_writeTABEAM
import test_dlpoly_writeTABLE

class TestPotentialFunctions(unittest.TestCase):
  """Tests for the potential function factories defined in atsim_potentials"""

  def testBuck(self):
    r = 0.5
    a = 1388.773
    rho = 2.76
    c = 175.0

    potfunc = atsim_potentials.buck(a,rho,c)
    self.assertAlmostEqual(-10041.34343, potfunc(r), places=5)

  def testHbnd(self):
    a = 300.0
    b = 20.0

    potfunc = atsim_potentials.hbnd(a,b)
    self.assertAlmostEqual(1208320, potfunc(0.5))


  def testPlus(self):
    potfunc = atsim_potentials.plus( atsim_potentials.buck(1388.773, 2.76, 175.0), atsim_potentials.hbnd(300.0, 20.0))
    self.assertAlmostEqual(1198278.65656831, potfunc(0.5))

class TestPotentialClass(unittest.TestCase):
  def testEnergy(self):
    """Check atsim_potentials.Potential.energy"""
    potfunc = atsim_potentials.buck(1388.773, 2.76, 175)
    pot = atsim_potentials.Potential("A", "B", potfunc)
    self.assertAlmostEquals(-10041.34343169, pot.energy(0.5), places = 5)

  def testForce(self):
    """Check atsim_potentials.Potential.force() method"""
    energyfunc = atsim_potentials.buck(1388.773, 0.362318841, 175.0)
    pot = atsim_potentials.Potential("A", "B", energyfunc)
    expect = 7.151
    actual = pot.force(2.0, h= 0.1e-5)
    self.assertAlmostEquals(expect, actual, places = 3)

class TestWritePotentials(unittest.TestCase):
  """Tests for the atsim_atsim_potentials.writePotentials"""

  def setUp(self):
    f_OO = atsim_potentials.potentialforms.buck(1633.00510, 0.327022, 3.948790)
    f_UU = atsim_potentials.potentialforms.buck(294.640000, 0.327022, 0.0)
    potential_objects = [
        atsim_potentials.Potential('O', 'O', f_OO),
        atsim_potentials.Potential('U', 'U', f_UU)]
    self.potential_objects = potential_objects

  def testDL_POLY(self):
    """Test DL_POLY tabulation"""
    sio = StringIO.StringIO()
    atsim_potentials.writePotentials(
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
    """Test atsim_potentials.writePotentials() for LAMMPS"""

    pota = atsim_potentials.Potential("A", "B", lambda x: x)
    potb = atsim_potentials.Potential("C", "D", lambda x: 6.0-x)

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
    atsim_potentials.writePotentials('LAMMPS', [pota,potb], 6.0, 6, sbuild)
    sbuild.seek(0)
    actual = sbuild.readlines()
    msg = "%s != %s" % (expect, actual)

    self.assertEquals(len(expect), len(actual), msg = msg)
    for e,a in zip(expect, actual):
      self.assertEquals(os.linesep, a[-1])
      a = a[:-1]
      self.assertEquals(e,a)
