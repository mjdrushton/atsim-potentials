import pytest
import py

from deepdiff import DeepDiff

from atsim.potentials._spline import SplinePotential
from atsim.potentials._spline import Exp_Spline
from atsim.potentials._spline import Spline_Point

from atsim.potentials import potentialforms
from atsim.potentials import plus
from atsim.potentials._util import gradient


# def _get_resource_dir(subdir):
#   rd = py.path.local(__file__).dirpath()
#   rd = rd.parts()[-1].join(subdir)
#   return rd

def test_exp_spline_spline_potential():
  bks_buck = potentialforms.buck(18003.7572, 0.2052048149, 133.5381)
  bks_coul = potentialforms.coul(2.4, -1.2)

  bks = plus(bks_buck, bks_coul)
  zbl = potentialforms.zbl(14, 8)
  spline = SplinePotential( zbl, bks_buck, 0.8, 1.4)

  expect_sc = (115.85848061468926, -554.8944421744643, 1103.3521708648423, -1086.8977963610287, 526.616672986912, -100.69329806728238, 0.0)
  actual = spline.splineCoefficients

  assert DeepDiff(expect_sc, actual) == {}

  for i in range(int(1.4 - 0.8 / 0.1)+1):
    r = float(i) * 0.1 + 0.8
    args = [r]
    args.extend(expect_sc)
    expect  = potentialforms.exp_spline(*args)
    assert pytest.approx(expect) == spline(r)

def test_exp_spline():
  bks_buck = potentialforms.buck(18003.7572, 1.0/4.87318, 133.5381)
  bks_coul = potentialforms.coul(2.4, -1.2)

  bks = plus(bks_buck, bks_coul)
  zbl = potentialforms.zbl(14, 8)
  
  r_detach = 0.8
  r_attach = 1.4

  v_detach = zbl(r_detach)
  v_attach = bks_buck(r_attach)

  dvdr_detach = gradient(zbl)
  dvdr_attach = gradient(bks_buck)

  ddvdr_detach = gradient(dvdr_detach)
  ddvdr_attach = gradient(dvdr_attach)

  detach = Spline_Point(zbl, r_detach)
  attach = Spline_Point(bks_buck, r_attach)

  assert pytest.approx(v_detach) == detach.v
  assert pytest.approx(v_attach) == attach.v

  assert pytest.approx(dvdr_detach(r_detach)) == detach.deriv
  assert pytest.approx(dvdr_attach(r_attach)) == attach.deriv

  assert pytest.approx(ddvdr_detach(r_detach)) == detach.deriv2
  assert pytest.approx(ddvdr_attach(r_attach)) == attach.deriv2

  spline = Exp_Spline(detach, attach)

  expect = (115.88493484897538, -555.0314373074318, 1103.6313912081537, -1087.1772275402761, 526.7537523874422, -100.71965947571076, 0.0)
  actual = spline.spline_coefficients
  assert DeepDiff(expect, actual) == {}




# def test_polynomial_spline():

#   left = potentialforms.bornmayer(11272.6, 0.1363)
#   right = 