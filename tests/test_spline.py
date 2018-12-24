import pytest
import py

from deepdiff import DeepDiff

from atsim.potentials.spline import SplinePotential
from atsim.potentials.spline import Exp_Spline
from atsim.potentials.spline import Spline_Point
from atsim.potentials.spline import Buck4_Spline

from atsim.potentials import potentialforms
from atsim.potentials import plus
from atsim.potentials._util import gradient


def test_exp_spline_spline_potential():
  bks_buck = potentialforms.buck(18003.7572, 0.2052048149, 133.5381)
  bks_coul = potentialforms.coul(2.4, -1.2)

  bks = plus(bks_buck, bks_coul)
  zbl = potentialforms.zbl(14, 8)
  spline = SplinePotential( zbl, bks_buck, 0.8, 1.4)

  expect_sc = (115.86040201825199, -554.9048594160315, 1103.3742860731818, -1086.9207282302543, 526.6282676867107, -100.6955851023172, 0.0)
  actual = spline.splineCoefficients

  assert DeepDiff(expect_sc, actual, significant_digits=5) == {}

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

  expect = (115.86040051844031, -554.9048516522953, 1103.3742702460838, -1086.9207123704682, 526.6282598808259, -100.69558359094829, 0.0)
  actual = spline.spline_coefficients
  assert DeepDiff(expect, actual, significant_digits=5) == {}

def test_buck4_spline():
  A = 11272.6
  rho = 0.1363
  C = 134.0

  born_mayer = potentialforms.bornmayer(A, rho)
  dispersion = potentialforms.buck(0, 1, C)

  r_detach = 1.2
  r_min = 2.1
  r_attach = 2.6

  detach = Spline_Point(born_mayer, r_detach)
  attach = Spline_Point(dispersion, r_attach)

  spline = Buck4_Spline(detach, attach, r_min)

  spline5_coeffs_expect = (479.9553, -1372.5304, 1562.2233, -881.9685,  246.4347,  -27.2447)
  spline3_coeffs_expect = (42.8917, -55.4965, 23.0774, -3.1314)

  assert pytest.approx(spline.spline5(1.3)) == 0.8227238
  DeepDiff(spline5_coeffs_expect, spline.spline5.args, significant_digits=4) == {}

  assert pytest.approx(spline.spline3(2.4)) == -0.662829
  DeepDiff(spline3_coeffs_expect, spline.spline3.args, significant_digits=4) == {}

  all_coefficients = spline5_coeffs_expect + spline3_coeffs_expect
  assert DeepDiff(all_coefficients, spline.spline_coefficients, significant_digits = 4) == {}

  assert pytest.approx(spline(1.0)) == 6.869673
  assert pytest.approx(spline(2.0)) == -0.8365282
  assert pytest.approx(spline(2.1)) == -0.8797423
  assert pytest.approx(spline(2.4)) == -0.662829
  assert pytest.approx(spline(2.6)) == -0.4337752
  assert pytest.approx(spline(2.8)) == -0.3125235
