"""Tests for the atsim.potentials.writePotentials"""

from io import StringIO

import pytest

from atsim import potentials

@pytest.fixture
def potential_objects():
  f_OO = potentials.potentialforms.buck(1633.00510, 0.327022, 3.948790)
  f_UU = potentials.potentialforms.buck(294.640000, 0.327022, 0.0)
  potential_objects = [
      potentials.Potential('O', 'O', f_OO),
      potentials.Potential('U', 'U', f_UU)]
  return potential_objects

def testDL_POLY(potential_objects):
  """Test DL_POLY tabulation"""
  sio = StringIO()
  potentials.writePotentials(
      'DL_POLY',
      potential_objects,
      6.5, 6500,
      out = sio)
  sio.seek(0)
  next(sio)
  line = next(sio)
  tokens = line.split()
  delpot, cutpot, ngrid = tokens

  assert pytest.approx(6.5/float(6500), abs = 1e-4) == float(delpot)
  assert pytest.approx(6.5) ==  float(cutpot)
  assert 6500 == int(ngrid)

  def testLAMMPS():
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

    sbuild = StringIO()
    potentials.writePotentials('LAMMPS', [pota,potb], 6.0, 6, sbuild)
    sbuild.seek(0)
    actual = sbuild.readlines()
    msg = "%s != %s" % (expect, actual)

    assert len(expect) == len(actual),  msg
    for e,a in zip(expect, actual):
      a = a.decode()
      assert os.linesep == a[-1]
      a = a[:-1]
      assert e == a
