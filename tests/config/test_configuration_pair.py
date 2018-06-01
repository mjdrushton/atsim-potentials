"""Tests for atsim.potential.config.Configuration for pair-potential models"""

import io
import os

import pytest

from atsim.potentials.config import Configuration
from atsim.potentials import potentialfunctions as pf

from .._runlammps import needsLAMMPS, extractLAMMPSEnergy, runLAMMPS, lammps_run_fixture
from .._rundlpoly import needsDLPOLY, runDLPoly, extractDLPOLYEnergy

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
  
[Pair]
O-O = as.buck 1000.0 0.3 32.0
U-O = as.buck 2000.0 0.2 0.0
"""

  cfgobj = Configuration()
  tabulation = cfgobj.read(io.StringIO(cfg_string))
  
  assert tabulation.nr == 1001
  assert tabulation.dr == 0.01
  assert tabulation.cutoff == 10.0

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