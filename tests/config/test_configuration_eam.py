import py
import pytest

import io

from atsim.potentials.config import Configuration
from atsim.potentials.config._config_parser import ConfigParser


from .._runlammps import needsLAMMPS, extractLAMMPSEnergy, runLAMMPS, lammps_run_fluorite_fixture
from .._rundlpoly import needsDLPOLY, runDLPoly, extractDLPOLYEnergy

def _get_resource_dir(subdir):
  rd = py.path.local(__file__).dirpath()
  rd = rd.parts()[-2].join(subdir)
  return rd

def _get_lammps_resource_dir():
  return _get_resource_dir('lammps_resources')

def _get_dlpoly_resource_dir():
  return _get_resource_dir('dl_poly_resources')


def test_configuration_setfl_synonyms():
  cfg_file_path = _get_lammps_resource_dir().join("CRG_U_Th.aspot")
  with cfg_file_path.open() as config_file:
    config_parser = ConfigParser(config_file)
  assert config_parser.tabulation.target == "setfl"

  # Now change setfl to lammps_eam_alloy and check that the target still registers as setfl
  from backports import configparser
  inifile = configparser.ConfigParser()
  inifile.read(cfg_file_path.strpath)

  inifile["Tabulation"]['target'] = "lammps_eam_alloy"

  modified = io.StringIO()
  inifile.write(modified)
  modified.seek(0)

  inifile = configparser.ConfigParser()
  inifile.read_file(modified)
  assert inifile["Tabulation"]['target'] == "lammps_eam_alloy"

  modified.seek(0)
  config_parser = ConfigParser(modified)
  assert config_parser.tabulation.target == "setfl"

@needsLAMMPS
def test_lammps_setfl_crg_tabulate_ThO2(lammps_run_fluorite_fixture):
  tmpdir = lammps_run_fluorite_fixture

  cfgobj = Configuration()
  config_file = _get_lammps_resource_dir().join("CRG_U_Th.aspot").open()
  tabulation = cfgobj.read(config_file)

  with tmpdir.join("table.eam.alloy").open("w") as outfile:
    tabulation.write(outfile)

  with lammps_run_fluorite_fixture.join("potentials.lmpinc").open('w') as potfile:
    potfile.write("""variable O equal 1
set type 1 charge -1.1104
set type 2 charge 2.2208

kspace_style pppm 1.0e-6

pair_style hybrid/overlay coul/long 10.0 eam/alloy
pair_coeff * * coul/long
pair_coeff * * eam/alloy table.eam.alloy O Th
""")

  runLAMMPS(cwd = tmpdir.strpath)
  energy = extractLAMMPSEnergy(cwd = tmpdir.strpath)

  expect = -157.552359260862
  assert pytest.approx(expect, rel = 1e-3) == energy  

@needsLAMMPS
def test_lammps_setfl_crg_tabulate_UO2(lammps_run_fluorite_fixture):
  tmpdir = lammps_run_fluorite_fixture

  cfgobj = Configuration()
  config_file = _get_lammps_resource_dir().join("CRG_U_Th.aspot").open()
  tabulation = cfgobj.read(config_file)

  with tmpdir.join("table.eam.alloy").open("w") as outfile:
    tabulation.write(outfile)

  with lammps_run_fluorite_fixture.join("potentials.lmpinc").open('w') as potfile:
    potfile.write("""variable O equal 1
set type 1 charge -1.1104
set type 2 charge 2.2208

kspace_style pppm 1.0e-6

pair_style hybrid/overlay coul/long 10.0 eam/alloy
pair_coeff * * coul/long
pair_coeff * * eam/alloy table.eam.alloy O U
""")

  runLAMMPS(cwd = tmpdir.strpath)
  energy = extractLAMMPSEnergy(cwd = tmpdir.strpath)

  expect = -162.708748563403
  assert pytest.approx(expect) == energy  

@needsDLPOLY
def test_dlpoly_TABEAM_tabulate_CeO2(tmpdir):
  # Copy files into the tmpdir.
  rd = _get_dlpoly_resource_dir()
  files = [
    ("CONFIG_CeO2", "CONFIG"),
    ("CONTROL_CeO2", "CONTROL"),
    ("FIELD_CeO2", "FIELD")
  ]

  for src, dest in files:
    src = rd.join(src)
    dest = tmpdir.join(dest)
    src.copy(dest)

  # Tabulate the TABEAM potential
  cfgobj = Configuration()
  config_file = rd.join("CRG_Ce.aspot").open()
  tabulation = cfgobj.read(config_file)

  with tmpdir.join("TABEAM").open("w") as outfile:
    tabulation.write(outfile)

  runDLPoly(cwd = tmpdir.strpath)
  actual = extractDLPOLYEnergy(cwd = tmpdir.strpath)

  expect = -532.6778
  assert pytest.approx(expect) == actual