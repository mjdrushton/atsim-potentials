"""Tests for atsim.potential.config.Configuration"""

import io
import os

import pytest

from atsim.potentials.config import Configuration
from atsim.potentials import potentialfunctions as pf

from .._runlammps import needsLAMMPS, extractLAMMPSEnergy, runLAMMPS, lammps_run_fixture

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
def test_pair_configuration_tabulate(lammps_run_fixture):
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
  energy = extractLAMMPSEnergy()
    

