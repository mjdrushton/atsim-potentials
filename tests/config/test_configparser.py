import io
import os

import pytest

from atsim.potentials.config import ConfigParser
from atsim.potentials.config import ConfigParserException
from atsim.potentials.config._config_parser import _RawConfigParser

from ._common import _get_lammps_resource_dir, _get_dlpoly_resource_dir


def test_potential_forms():
  """Testing reading of potential forms from [Potential-Form]"""

  cfg_string = u"""[Potential-Form]
buck_morse(r_ij, A,rho,C,D,gamma,r0) : buck(r_ij, A,rho,C) + morse(r_ij, gamma,r0,D)
density(r_ij, n) : (n/r_ij^8) * (1/2)*(1+erf(20*(r_ij-1.5)))
"""

  expect = [
    (("buck_morse", ["r_ij", "A", "rho", "C", "D", "gamma", "r0"]), "buck(r_ij, A,rho,C) + morse(r_ij, gamma,r0,D)"),
    (("density", ["r_ij", "n"]), "(n/r_ij^8) * (1/2)*(1+erf(20*(r_ij-1.5)))")]

  parsed = ConfigParser(io.StringIO(cfg_string))
  actual = parsed.potential_form
  assert expect == actual

def test_parse_potential_form_signature():
  pfstr = "buck_morse(r_ij, A,rho,C,D,gamma,r0)"

  cfg = ConfigParser(io.StringIO())
  expect = ("buck_morse", ["r_ij", "A", "rho", "C", "D", "gamma", "r0"])
  actual = cfg._parse_potential_form_signature(pfstr)
  assert expect == actual

def test_empty_tabulation_section():
  cfg_string = u""
  parsed = ConfigParser(io.StringIO(cfg_string))
  assert not parsed.tabulation is None


def test_tabulation_cutoff_and_step():

  cfg_string = u"""[Tabulation]
cutoff : 10.0
"""

  parsed = ConfigParser(io.StringIO(cfg_string))
  assert 10.0 == parsed.tabulation.cutoff
  assert parsed.tabulation.nr is None
  
  cfg_string = u"""[Tabulation]
cutoff : 10.0
nr : 1000
"""

  parsed = ConfigParser(io.StringIO(cfg_string))
  assert 10.0 == parsed.tabulation.cutoff
  assert parsed.tabulation.nr == 1000

  cfg_string = u"""[Tabulation]
cutoff : 10.0
dr : 1.0
"""

  parsed = ConfigParser(io.StringIO(cfg_string))
  assert 10.0 == parsed.tabulation.cutoff
  assert parsed.tabulation.nr == 11

  cfg_string = u"""[Tabulation]
nr : 11
dr : 1.0
"""

  parsed = ConfigParser(io.StringIO(cfg_string))
  assert 10.0 == parsed.tabulation.cutoff
  assert parsed.tabulation.nr == 11

  with pytest.raises(ConfigParserException):
    cfg_string = u"""[Tabulation]
dr : 11
"""
    parsed = ConfigParser(io.StringIO(cfg_string))
    parsed.tabulation
    
  with pytest.raises(ConfigParserException):
    cfg_string = u"""[Tabulation]
cutoff : 10.0
dr : 0.001
nr : 11
"""
    parsed = ConfigParser(io.StringIO(cfg_string))
    parsed.tabulation

def test_tabulation_density_cutoff_and_step():

  cfg_string = u"""[Tabulation]
cutoff_rho : 10.0
"""

  parsed = ConfigParser(io.StringIO(cfg_string))
  assert 10.0 == parsed.tabulation.cutoff_rho
  assert parsed.tabulation.nrho is None
  
  cfg_string = u"""[Tabulation]
cutoff_rho : 10.0
nrho : 1000
"""

  parsed = ConfigParser(io.StringIO(cfg_string))
  assert 10.0 == parsed.tabulation.cutoff_rho
  assert parsed.tabulation.nrho == 1000

  cfg_string = u"""[Tabulation]
cutoff_rho : 10.0
drho : 1.0
"""

  parsed = ConfigParser(io.StringIO(cfg_string))
  assert 10.0 == parsed.tabulation.cutoff_rho
  assert parsed.tabulation.nrho == 11

  cfg_string = u"""[Tabulation]
nrho : 11
drho : 1.0
"""

  parsed = ConfigParser(io.StringIO(cfg_string))
  assert 10.0 == parsed.tabulation.cutoff_rho
  assert parsed.tabulation.nrho == 11

  with pytest.raises(ConfigParserException):
    cfg_string = u"""[Tabulation]
drho : 11
"""
    parsed = ConfigParser(io.StringIO(cfg_string))
    parsed.tabulation
    
  with pytest.raises(ConfigParserException):
    cfg_string = u"""[Tabulation]
cutoff_rho : 10.0
drho : 0.001
nrho : 11
"""
    parsed = ConfigParser(io.StringIO(cfg_string))
    parsed.tabulation

def test_parsed_sections():
  expect = ['tabulation', 'potential_form', 'eam_embed', 'eam_density', 'pair']
  expect.sort()

  with _get_dlpoly_resource_dir().join("CRG_Ce.aspot").open() as infile:
    cp = ConfigParser(infile)
    actual = cp.parsed_sections
    actual.sort()
  assert expect == actual

  expect = ['tabulation', 'potential_form', 'pair']
  expect.sort()

  with _get_lammps_resource_dir().join("zbl_spline.aspot").open() as infile:
    cp = ConfigParser(infile)
    actual = cp.parsed_sections
    actual.sort()
  assert expect == actual

  expect = ['tabulation', 'potential_form', 'eam_embed', 'eam_density_fs', 'pair']
  expect.sort()

  with _get_lammps_resource_dir().join("AlFe_setfl_fs.aspot").open() as infile:
    cp = ConfigParser(infile)
    actual = cp.parsed_sections
    actual.sort()
  assert expect == actual

def test_orphan_sections():
  with _get_lammps_resource_dir().join("AlFe_setfl_fs.aspot").open("r") as infile:
    sio = io.StringIO(infile.read())

  sio.seek(0,os.SEEK_END)
  sio.write("\n\n[Charges]\nU : 2.22\nO: -1.11\n")
  sio.write("[Masses]\nU : 1.234\nO: 4.56\n")
  sio.seek(0)

  expect = ['tabulation', 'potential_form', 'eam_embed', 'eam_density_fs', 'pair']
  expect.sort()

  cp = ConfigParser(sio)
  actual = cp.parsed_sections
  actual.sort()
  assert expect == actual

  actual = cp.orphan_sections
  expect = ["Charges", "Masses"]

  actual.sort()
  expect.sort()

  assert expect == actual

  



