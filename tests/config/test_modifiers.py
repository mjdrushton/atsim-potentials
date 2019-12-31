import io

import pytest

from atsim.potentials import potentialfunctions as pf
from atsim.potentials.config import Configuration
from atsim.potentials.config._common import Modifier_Exception, Unknown_Modifier_Exception

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


def test_pow_modifier():
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

def test_pow_modifier():
  cfg_string = u"""[Pair]

A-B = pow(as.polynomial 3.0 2.0, as.constant 2)
A-A = pow(as.polynomial 3.0 2.0, pow(as.polynomial 0 2, as.constant 0.5))
A-C = pow(sum(as.buck 1000 0.3 0, as.polynomial 0 1 2), pow(sum(as.polynomial 0 0.1, as.polynomial 0 -0.05), as.constant 0.5))
A-D = pow(as.polynomial 1.0 2.0 3.0, as.polynomial 0 5 0.1, as.constant 0.01)

"""

  cfgobj = Configuration()
  tabulation = cfgobj.read(io.StringIO(cfg_string))

  pots = tabulation.potentials
  potdict = dict([((p.speciesA, p.speciesB),p) for p in pots])
  ab = potdict[("A", "B")]

  import sympy
  r,A,B,C = sympy.symbols("r A B C")

  ab_sympy = (A + B*r)**C
  var_dict = dict(r = 2.5, A = 3.0, B = 2.0, C = 2.0)

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


  ac = potdict[("A", "A")]
  r = sympy.symbols("r")
  # A-C = pow(sum(as.buck 100 0.2 0, as.polynomial 0 1 2), pow(sum(as.polynomial 0 0.1, as.polynomial 0 -0.05), as.constant 0.5))
  # pow(as.polynomial 3.0 2.0, pow(as.constant 3, as.constant 2))
  ac_sympy = (3.0 + 2.0*r)**((2*r)**0.5)
  var_dict = dict(r=2.5)
  actual = ac.energy(var_dict["r"])
  expect = float(ac_sympy.subs(var_dict))
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


  ac = potdict[("A", "C")]

  r = sympy.symbols("r")
  A_buck, rho_buck, C_buck = sympy.symbols("A_buck rho_buck C_buck")

  # A-C = pow(sum(as.buck 100 0.2 0, as.polynomial 0 1 2), pow(sum(as.polynomial 0 0.1, as.polynomial 0 -0.05), as.constant 0.5))
  ac_sympy = ((A_buck * sympy.exp(- r/rho_buck)) + (r + 2*r**2))**((0.1*r -0.05*r)**(0.5))
  var_dict = dict(r=2.5, A_buck = 1000, rho_buck = 0.3)

  actual = ac.energy(var_dict["r"])
  expect = float(ac_sympy.subs(var_dict))
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

  ac = potdict[("A", "D")]

  # A-D = pow(as.polynomial 1.0 2.0 3.0, as.polynomial 4 5 6, as.constant 3)

  r = sympy.symbols("r")
  p1, p2, p3, p4, p5, p6, p7 = sympy.symbols("p1 p2 p3 p4 p5 p6 p7")

  ac_sympy = ((p1 + p2*r + p3*r**2)**(p4 + p5*r + p6*r**2))**p7
  var_dict = dict(r=2.5, 
  p1 = 1, p2 = 2, p3 = 3, p4 = 0, p5 = 5, p6 = 0.1, p7 = 1.0/100.0)

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

def test_trans_modifier():
  cfg_string = u"""[Pair]

A-B = trans(pow(as.polynomial 3.0 2.0, as.constant 2), as.constant -2.0)
A-A = sum(as.buck 1000.0 0.1 0, trans(as.buck 1000.0 0.1 0, as.constant 1.0))
"""

  cfgobj = Configuration()
  tabulation = cfgobj.read(io.StringIO(cfg_string))

  pots = tabulation.potentials
  potdict = dict([((p.speciesA, p.speciesB),p) for p in pots])
  ab = potdict[("A", "B")]

  import sympy
  r,A,B,C = sympy.symbols("r A B C")

  ab_sympy = (A + B*(r-2.0))**C
  var_dict = dict(r = 2.5, A = 3.0, B = 2.0, C = 2.0)

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

  aa = potdict[("A", "A")]
  r,A,B,C,D,E,F = sympy.symbols("r A B C D E F")
  aa_sympy = ( A * sympy.exp(-r/B) - C/r**6 ) + ( D * sympy.exp(-(r+1)/E) - F/(r+1)**6 )
  var_dict = dict(r = 2.5, A = 1000.0, B = 0.1, C = 0.0, D=1000.0, E=0.1, F=0.0)

  expect = float(aa_sympy.subs(var_dict))
  actual = aa.energy(var_dict["r"])
  assert pytest.approx(expect) == actual

  assert hasattr(aa.potentialFunction, "deriv")
  aa_deriv = sympy.diff(aa_sympy, r)
  expect = float(aa_deriv.subs(var_dict))
  actual = aa.potentialFunction.deriv(var_dict["r"])
  assert pytest.approx(expect) == actual

  assert hasattr(aa.potentialFunction, "deriv2")
  aa_deriv2 = sympy.diff(aa_sympy, r,2)
  expect = float(aa_deriv2.subs(var_dict))
  actual = aa.potentialFunction.deriv2(var_dict["r"])
  assert pytest.approx(expect) == actual



def test_unknown_modifier_exception():
  cfg_string = u"""[Pair]
U-O = unknown_modifier(as.buck 1000.0 0.3 32.0, as.constant 1.0)
"""

  cfgobj = Configuration()
  with pytest.raises(Unknown_Modifier_Exception):
    tabulation = cfgobj.read(io.StringIO(cfg_string))
  