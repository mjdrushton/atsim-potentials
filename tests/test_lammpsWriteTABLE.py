"""Module containing tests for the _lammps_writeTABLE module"""

import unittest
import StringIO
import os

from atsim import potentials
from atsim.potentials import _lammps_writeTABLE

class LammpsWriteTABLETestCase(unittest.TestCase):
  """Test case for lammps.writeTABLE module"""

  def testWriteSinglePotential(self):
    """Test that atsim.potentials._lammps_writeTABLE._writeSinglePotential() works correctly"""

    expect=["A-B",
            "N 6 R 0.10000000 5.10000000",
            "",
            "1 0.10000000 4.90000000 1.00000000",
            "2 1.10000000 3.90000000 1.00000000",
            "3 2.10000000 2.90000000 1.00000000",
            "4 3.10000000 1.90000000 1.00000000",
            "5 4.10000000 0.90000000 1.00000000",
            "6 5.10000000 -0.10000000 1.00000000"]

    energyfunc = lambda x: 5.0-x
    pot =potentials.Potential("A", "B", energyfunc)

    sbuild = StringIO.StringIO()
    _lammps_writeTABLE._writeSinglePotential(pot, 0.1, 5.1, 6, sbuild)
    sbuild.seek(0)
    actual = sbuild.readlines()

    msg = "%s != %s" % (expect, actual)

    self.assertEquals(len(expect), len(actual), msg = msg)
    for e,a in zip(expect, actual):
      self.assertEquals(os.linesep, a[-1])
      a = a[:-1]
      self.assertEquals(e,a)

  def testWritePotentials(self):
    """Test _lammps_writeTABLE.writePotentials() function"""

    pota =potentials.Potential("A", "B", lambda x: x)
    potb =potentials.Potential("C", "D", lambda x: 5.0-x)

    expect=["A-B",
            "N 6 R 0.10000000 5.10000000",
            "",
            "1 0.10000000 0.10000000 -1.00000000",
            "2 1.10000000 1.10000000 -1.00000000",
            "3 2.10000000 2.10000000 -1.00000000",
            "4 3.10000000 3.10000000 -1.00000000",
            "5 4.10000000 4.10000000 -1.00000000",
            "6 5.10000000 5.10000000 -1.00000000",
            "",
            "C-D",
            "N 6 R 0.10000000 5.10000000",
            "",
            "1 0.10000000 4.90000000 1.00000000",
            "2 1.10000000 3.90000000 1.00000000",
            "3 2.10000000 2.90000000 1.00000000",
            "4 3.10000000 1.90000000 1.00000000",
            "5 4.10000000 0.90000000 1.00000000",
            "6 5.10000000 -0.10000000 1.00000000"]

    sbuild = StringIO.StringIO()
    _lammps_writeTABLE.writePotentials([pota,potb], 0.1, 5.1, 6, sbuild)
    sbuild.seek(0)
    actual = sbuild.readlines()

    msg = "%s != %s" % (expect, actual)

    self.assertEquals(len(expect), len(actual), msg = msg)
    for e,a in zip(expect, actual):
      self.assertEquals(os.linesep, a[-1])
      a = a[:-1]
      self.assertEquals(e,a)
