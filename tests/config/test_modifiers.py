import io

import pytest

from atsim.potentials import potentialfunctions as pf
from atsim.potentials.config import Configuration


def test_sum_modifier():
  cfg_string = u"""[Pair]
O-O = as.buck 1000.0 0.3 32.0
U-O = sum(as.buck 1000.0 0.3 32.0, as.constant 1.0)
"""

  cfgobj = Configuration()
  tabulation = cfgobj.read(io.StringIO(cfg_string))
  
  potlist = tabulation.potentials

  expect = [ ("O", "O"), ("U", "O")]
  actual = [(p.speciesA, p.speciesB) for p in potlist]
  assert sorted(expect) == sorted(actual)

  r = 1.3
  buck_oo = pf.buck(r, 1000.0, 0.3, 32.0)

  expect = [ 
    (("O", "O"), pytest.approx(buck_oo)), 
    (("U", "O"), pytest.approx(buck_oo + 1.0))]
  actual = [((p.speciesA, p.speciesB),p.energy(r)) for p in potlist]
  assert sorted(expect) == sorted(actual)

def test_sum_modifier_deriv():
  """Make sure the sum() modifier uses .deriv() methods correctly if they are available"""
  
  cfg_string = u"""[Pair]
#no deriv
A-B = sum(born_mayer 1000.0 0.1, dispersion 32.0)
# both have deriv
B-C = sum(as.bornmayer 1000.0 0.1,  as.buck 0 1.0 32.0)
#one_deriv
C-D = sum(as.bornmayer 1000.0 0.1, dispersion 32.0)

[Potential-Form]
born_mayer(r, A, rho) = A * exp(-r/rho) 
dispersion(r, C) = - C/r^6
  
"""
  cfgobj = Configuration()
  tabulation = cfgobj.read(io.StringIO(cfg_string))
  
  pots = tabulation.potentials

  expect = [ ("A", "B"), ("B", "C"), ("C", "D")]
  actual = [(p.speciesA, p.speciesB) for p in pots]
  assert sorted(expect) == sorted(actual)

  potdict = dict([((p.speciesA, p.speciesB),p) for p in pots])
  assert not hasattr(potdict[('A', 'B')].potentialFunction, 'deriv')
  assert  hasattr(potdict[('B', 'C')].potentialFunction, 'deriv')
  assert  hasattr(potdict[('C', 'D')].potentialFunction, 'deriv')

  expect_energy = pytest.approx(pf.buck(1.6, 1000, 0.1, 32.0))
  expect_deriv = pytest.approx(pf.buck.deriv(1.6, 1000, 0.1, 32.0))

  assert expect_energy == potdict[('A', 'B')].potentialFunction(1.6)
  
  assert expect_energy == potdict[('B', 'C')].potentialFunction(1.6)
  assert expect_deriv == potdict[('B', 'C')].potentialFunction.deriv(1.6)

  assert expect_energy == potdict[('C', 'D')].potentialFunction(1.6)
  assert expect_deriv == potdict[('C', 'D')].potentialFunction.deriv(1.6)

def test_product_modifier():
  cfg_string = u"""[Pair]

A-B = product(as.buck 1000.0 0.2 32.0, as.polynomial 0.0 2.0, as.polynomial 1.0 -3.0 0.5)
A-C = sum(
          product(
            as.constant 2.0, 
            as.polynomial 1.0 2.0 1.5, 
            product(
              as.polynomial 1.0 2.0 3.0, 
              as.polynomial 4.0 5.0 6.0
            ) ), 
          as.polynomial 0.0 -1.0)

"""

  cfgobj = Configuration()
  tabulation = cfgobj.read(io.StringIO(cfg_string))

  pots = tabulation.potentials
  potdict = dict([((p.speciesA, p.speciesB),p) for p in pots])
  ab = potdict[("A", "B")]

  import sympy
  r,A,rho,C,p1_0,p1_1,p2_0,p2_1,p2_2 = sympy.symbols("r A rho C p1_0 p1_1 p2_0 p2_1 p2_2")

  ab_sympy = (A * sympy.exp(-r/rho) - C/r**6) * (p1_0 + p1_1*r) * (p2_0 + p2_1*r + p2_2*r**2)
  var_dict = dict(r = 2.5, 
                  A = 1000.0, rho = 0.2, C = 32.0,
                  p1_0 = 0.0, p1_1  = 2.0,
                  p2_0 = 1.0, p2_1 = -3.0, p2_2 = 0.5)

  expect = float(ab_sympy.subs(var_dict))
  actual = ab.energy(var_dict["r"])
  assert pytest.approx(expect) == actual

  assert hasattr(ab.potentialFunction, "deriv")
  ab_deriv = sympy.diff(ab_sympy, r)
  expect = float(ab_deriv.subs(var_dict))
  actual = ab.potentialFunction.deriv(var_dict["r"])
  assert pytest.approx(expect) == actual

  assert hasattr(ab.potentialFunction, "deriv2")
  ab_deriv2 = sympy.diff(ab_sympy, r,2)
  expect = float(ab_deriv2.subs(var_dict))
  actual = ab.potentialFunction.deriv2(var_dict["r"])
  assert pytest.approx(expect) == actual

  ac = potdict[("A", "C")]

  r = sympy.symbols("r")
  A = sympy.symbols("A")
  p1_0, p1_1, p1_2 = sympy.symbols("p1_0 p1_1 p1_2")
  p2_0, p2_1, p2_2 = sympy.symbols("p2_0 p2_1 p2_2")
  p3_0, p3_1, p3_2 = sympy.symbols("p3_0 p3_1 p3_2")
  p4_0, p4_1 = sympy.symbols("p4_0 p4_1")

  ac_sympy = (A * (p1_0 + p1_1*r + p1_2*r**2) * ((p2_0 + p2_1 *r + p2_2*r**2) * (p3_0 + p3_1*r + p3_2*r**2))) + (p4_0 + p4_1*r)
  var_dict = dict(r=2.5, 
                  A = 2.0, 
                  p1_0 = 1.0, p1_1 = 2.0, p1_2 = 1.5,
                  p2_0 = 1.0, p2_1 = 2.0, p2_2 = 3.0,
                  p3_0 = 4.0, p3_1 = 5.0, p3_2 = 6.0,
                  p4_0 = 0.0, p4_1 = -1.0)

  expect = float(ac_sympy.subs(var_dict))
  actual = ac.energy(var_dict["r"])
  assert pytest.approx(expect) == actual

  assert hasattr(ac.potentialFunction, "deriv")
  ac_deriv = sympy.diff(ac_sympy, r)
  expect = float(ac_deriv.subs(var_dict))
  actual = ac.potentialFunction.deriv(var_dict["r"])
  assert pytest.approx(expect) == actual

  assert hasattr(ac.potentialFunction, "deriv2")
  ac_deriv2 = sympy.diff(ac_sympy, r,2)
  expect = float(ac_deriv2.subs(var_dict))
  actual = ac.potentialFunction.deriv2(var_dict["r"])
  assert pytest.approx(expect) == actual