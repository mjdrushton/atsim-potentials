"""Tests for atsim.potential.config.Configuration for pair-potential models"""

import io
import os

import pytest

from atsim.potentials.config import Configuration
from atsim.potentials import potentialfunctions as pf

from .._runlammps import needsLAMMPS, extractLAMMPSEnergy, runLAMMPS, lammps_run_fixture
from .._rundlpoly import needsDLPOLY, runDLPoly, extractDLPOLYEnergy
from .._rungulp import needsGULP, runGULP, extractGULPEnergy

def test_pair_configuration():
  cfg_string = u"""[Pair]
O-O = as.buck 1000.0 0.3 32.0
U-O = as.buck 2000.0 0.2 0.0
"""

  cfgobj = Configuration()
  tabulation = cfgobj.read(io.StringIO(cfg_string))
  
  assert tabulation.nr == 1001
  assert tabulation.dr == 0.01
  assert tabulation.cutoff == 10.0

  assert tabulation.type == "Pair"
  assert tabulation.target == "LAMMPS"

  assert 2 == len(tabulation.potentials)
  
  potlist = tabulation.potentials

  expect = [ ("O", "O"), ("U", "O")]
  actual = [(p.speciesA, p.speciesB) for p in potlist]
  assert sorted(expect) == sorted(actual)

  r = 1.3
  buck_oo = pf.buck(r, 1000.0, 0.3, 32.0)
  buck_uo = pf.buck(r, 2000.0, 0.2, 0.0)

  expect = [ 
    (("O", "O"), pytest.approx(buck_oo)), 
    (("U", "O"), pytest.approx(buck_uo))]
  actual = [((p.speciesA, p.speciesB),p.energy(r)) for p in potlist]
  assert sorted(expect) == sorted(actual)

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


def test_multirange_deriv():
  """Make sure that .deriv() methods are used where possible in MultiRangePotential"""
  # Neither potentialform has deriv() - resultant doesn't have deriv()
  cfg_string = u"""[Pair]
#no deriv
A-B = born_mayer 1000.0 0.1 >=3.0 dispersion 32.0
# both have deriv
B-C = as.bornmayer 1000.0 0.1 >=3.0 as.buck 0 1.0 32.0
#one_deriv
C-D = as.bornmayer 1000.0 0.1 >=3.0 dispersion 32.0

[Potential-Form]
born_mayer(r, A, rho) = A * exp(-r/rho) 
dispersion(r, C) = - C/r^6
  
"""

  cfgobj = Configuration()

  tabulation = cfgobj.read(io.StringIO(cfg_string))
  pots = tabulation.potentials

  potdict = dict([((p.speciesA, p.speciesB),p) for p in pots])

  assert not hasattr(potdict[('A', 'B')].potentialFunction, 'deriv')
  assert  hasattr(potdict[('B', 'C')].potentialFunction, 'deriv')
  assert  hasattr(potdict[('C', 'D')].potentialFunction, 'deriv')

  assert pytest.approx(pf.bornmayer.deriv(1.6, 1000.0, 0.1)) == potdict[('B', 'C')].potentialFunction.deriv(1.6)
  assert pytest.approx(pf.buck.deriv(3.2, 0, 1.0, 32.0)) == potdict[('B', 'C')].potentialFunction.deriv(3.2)

  assert pytest.approx(pf.bornmayer.deriv(1.6, 1000.0, 0.1)) == potdict[('C', 'D')].potentialFunction.deriv(1.6)
  assert pytest.approx(pf.buck.deriv(3.2, 0, 1.0, 32.0)) == potdict[('C', 'D')].potentialFunction.deriv(3.2)

def test_pair_deriv():
  """Test that pair potentials defined in potentialfunctions expose their .deriv() functions"""

  cfg_string = u"""[Pair]
O-O = as.buck 1000.0 0.3 32.0
"""

  cfgobj = Configuration()

  # import pdb; pdb.set_trace()
  tabulation = cfgobj.read(io.StringIO(cfg_string))

  potfunc = tabulation.potentials[0].potentialFunction
  assert hasattr(potfunc, "deriv")
  assert hasattr(potfunc, "deriv2")

  assert pf.buck.deriv(1.6, 1000.0, 0.3, 32.0) == potfunc.deriv(1.6)
  assert pf.buck.deriv2(1.6, 1000.0, 0.3, 32.0) == potfunc.deriv2(1.6)

  # Test that cexprtk potentialform doesn't provide a .deriv function.
  cfg_string = u"""[Pair]
O-O = my_buck 1000.0 0.3 32.0

[Potential-Form]
my_buck(r, A, rho, C) = A*exp(-r/rho) - C/r^6

"""

  cfgobj = Configuration()

  # import pdb; pdb.set_trace()
  tabulation = cfgobj.read(io.StringIO(cfg_string))

  potfunc = tabulation.potentials[0].potentialFunction
  assert not hasattr(potfunc, "deriv")
  assert not hasattr(potfunc, "deriv2")

def test_vararg_potentialform():
  """Check that a potential form that uses varargs works correctly"""

  cfg_string = u"""[Pair]
O-O = as.polynomial 10.0 20.0 30.0
"""

  cfgobj = Configuration()
  tabulation = cfgobj.read(io.StringIO(cfg_string))

  potfunc = tabulation.potentials[0].potentialFunction
  assert pf.polynomial(1.6, 10.0, 20.0, 30.0) == potfunc(1.6)

  # Now check its use in a custom potentialform
  cfg_string = u"""[Pair]
O-O = mypoly 10.0 20.0 30.0

[Potential-Form]
mypoly(r, A, B, C) = as.polynomial(r, A,B,C) + 10.0"""

  expect = pf.polynomial(1.6, 10.0, 20.0, 30.0) + 10.0

  cfgobj = Configuration()
  tabulation = cfgobj.read(io.StringIO(cfg_string))

  potfunc = tabulation.potentials[0].potentialFunction
  assert pytest.approx(expect) == potfunc(1.6)

@needsLAMMPS
def test_lammps_pair_configuration_tabulate(lammps_run_fixture):
  tmpdir = lammps_run_fixture
  cfg_string = u"""[Potential-Form]
buck(r, A, rho, C) : A*exp(-r/rho) - C/r^6
morse(r, gamma, r_star, D) : D*(exp(-2.0*gamma*(r-r_star)) - 2.0*exp(-gamma*(r-r_star)))
buck_morse(r, A, rho, C, gamma, r_star, D) : buck(r,A,rho,C) + morse(r, gamma, r_star, D)

[Pair]
O-U = buck_morse 1000.0 0.1 32.0 0.3 2.0 10.0

"""

  cfgobj = Configuration()
  tabulation = cfgobj.read(io.StringIO(cfg_string))

  with tmpdir.join("table.lmptab").open("w") as outfile:
    tabulation.write(outfile)

  with lammps_run_fixture.join("potentials.lmpinc").open('w') as potfile:
    potfile.write("pair_style hybrid/overlay zero {0} table linear {1}\n".format(tabulation.cutoff, tabulation.nr))
    potfile.write("pair_coeff * * zero\n")
    potfile.write("pair_coeff 1 2 table table.lmptab O-U\n")
    potfile.write("\n")

  runLAMMPS(cwd = tmpdir.strpath)
  energy = extractLAMMPSEnergy(cwd = tmpdir.strpath)

  expect = pf.buck(2.0, 1000.0, 0.1, 32.0) + pf.morse(2.0, 0.3, 2.0, 10.0)
  assert pytest.approx(expect) == energy

def test_dlpoly_pair_configuration():
  cfg_string = u"""[Tabulation]
target : DLPOLY
nr : 1000
cutoff : 9.99
  
[Pair]
O-O = as.buck 1000.0 0.3 32.0
U-O = as.buck 2000.0 0.2 0.0
"""

  cfgobj = Configuration()
  tabulation = cfgobj.read(io.StringIO(cfg_string))
  
  assert tabulation.nr == 1000
  assert tabulation.dr == 0.01
  assert tabulation.cutoff == 9.99

  assert tabulation.type == "Pair"
  assert tabulation.target == "DLPOLY"

  assert 2 == len(tabulation.potentials)
  
  potlist = tabulation.potentials

  expect = [ ("O", "O"), ("U", "O")]
  actual = [(p.speciesA, p.speciesB) for p in potlist]
  assert sorted(expect) == sorted(actual)

  r = 1.3
  buck_oo = pf.buck(r, 1000.0, 0.3, 32.0)
  buck_uo = pf.buck(r, 2000.0, 0.2, 0.0)

  expect = [ 
    (("O", "O"), pytest.approx(buck_oo)), 
    (("U", "O"), pytest.approx(buck_uo))]
  actual = [((p.speciesA, p.speciesB),p.energy(r)) for p in potlist]
  assert sorted(expect) == sorted(actual)

@needsDLPOLY
def test_dlpoly_pair_configuration_tabulate(tmpdir):
  cfg_string = u"""[Tabulation]
target : DLPOLY
nr : 1000
cutoff : 6.5

[Potential-Form]
buck(r, A, rho, C) : A*exp(-r/rho) - C/r^6
morse(r, gamma, r_star, D) : D*(exp(-2.0*gamma*(r-r_star)) - 2.0*exp(-gamma*(r-r_star)))
buck_morse(r, A, rho, C, gamma, r_star, D) : buck(r,A,rho,C) + morse(r, gamma, r_star, D)

[Pair]
O-U = buck_morse 1000.0 0.1 32.0 0.3 2.0 10.0

"""
  import py.path

  this_dir = py.path.local(py.path.local(__file__).dirname)
  dlpoly_resource_dir = py.path.local(this_dir.dirname).join("dl_poly_resources")

  templatevars = dict(
    speciesA = "O",
    speciesB = "U",
    Ax = 0.0,
    Ay = 0.0,
    Az = 0.0,
    Bx = 2.0,
    By = 0.0,
    Bz = 0.0,
    potDef = "vdw 1\nO U tab\n")

  # Write the config file
  with dlpoly_resource_dir.join("CONFIG_pair.in").open() as config_template:
    with tmpdir.join("CONFIG").open("w") as outfile:
      outfile.write(config_template.read() % templatevars)

  # Write the FIELD file
  with dlpoly_resource_dir.join("FIELD_pair.in").open() as field_template:
    with tmpdir.join("FIELD").open("w") as outfile:
      outfile.write(field_template.read() % templatevars)

  # Write the CONTROL file
  dlpoly_resource_dir.join("CONTROL_pair").copy(tmpdir.join("CONTROL"))

  # Finally write the TABLE file
  cfgobj = Configuration()
  tabulation = cfgobj.read(io.StringIO(cfg_string))

  with tmpdir.join("TABLE").open("w") as outfile:
    tabulation.write(outfile)

  runDLPoly(cwd = tmpdir.strpath)
  energy = extractDLPOLYEnergy(cwd = tmpdir.strpath)

  expect = pf.buck(2.0, 1000.0, 0.1, 32.0) + pf.morse(2.0, 0.3, 2.0, 10.0)
  assert pytest.approx(expect) == energy

@needsGULP
def test_gulp_pair_configuration_tabulate(tmpdir):
  cfg_string = u"""[Tabulation]
target : GULP
dr: 0.01
cutoff : 15.0

[Pair]
O-U = as.bornmayer 1761.775 0.35642
O-O = as.buck 9547.96 0.2192 32.0

"""

  gulp_input = u"""single

cell
5.468 5.468 5.468 90.0 90.0 90.0

frac
U 0 0 0
U 1/2 1/2 0
U 1/2 0 1/2
U 0 1/2 1/2

O 1/4 1/4 1/4
O 1/4 3/4 1/4
O 3/4 3/4 1/4
O 3/4 1/4 1/4

O 1/4 1/4 3/4
O 1/4 3/4 3/4
O 3/4 3/4 3/4
O 3/4 1/4 3/4


species
U 4.0
O -2.0

include potentials.lib

"""
  
  # First calculate the expected energy using GULP's built-in analytical potentials
  with tmpdir.join("potentials.lib").open("w") as potfile:
    potfile.write("buck\n")
    potfile.write("O O 9547.96 0.2192 32.0 15.0\n")
    potfile.write("O U 1761.775 0.35642 0.0 15.0\n")

  gulp_infile = io.StringIO(gulp_input)
  gulp_infile.seek(0)

  gulp_outfile = io.StringIO()
  runGULP(gulp_infile, gulp_outfile, cwd = tmpdir.strpath)

  gulp_outfile.seek(0)
  expect = extractGULPEnergy(gulp_outfile)

  tmpdir.join("potentials.lib").remove()
  assert not tmpdir.join("potentials.lib").exists()

  # Now build a potential model and tabulate it - then re-run the calculation and check the energies match.
  cfgobj = Configuration()
  tabulation = cfgobj.read(io.StringIO(cfg_string))

  with tmpdir.join("potentials.lib").open("w") as potfile:
    tabulation.write(potfile)

  gulp_infile.seek(0)

  gulp_outfile = io.StringIO()
  runGULP(gulp_infile, gulp_outfile, cwd = tmpdir.strpath)

  gulp_outfile.seek(0)
  actual = extractGULPEnergy(gulp_outfile)
  assert pytest.approx(expect, abs=1e-3) == actual
