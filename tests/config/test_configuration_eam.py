import py
import pytest

import io

from atsim.potentials.config import Configuration
from atsim.potentials.config._config_parser import ConfigParser


from .._runlammps import needsLAMMPS, extractLAMMPSEnergy, runLAMMPS, lammps_run_fluorite_fixture, lammps_run_fixture
from .._rundlpoly import needsDLPOLY, runDLPoly, extractDLPOLYEnergy

from ._common import _get_dlpoly_resource_dir, _get_lammps_resource_dir

def test_configuration_setfl_synonyms():
  cfg_file_path = _get_lammps_resource_dir().join("CRG_U_Th.aspot")
  with io.open(cfg_file_path.strpath, encoding = "utf8") as config_file:
    config_parser = ConfigParser(config_file)
  assert config_parser.tabulation.target == u"setfl"

  # Now change setfl to lammps_eam_alloy and check that the target still registers as setfl
  import configparser
  inifile = configparser.ConfigParser()
  inifile.read(u"{}".format(cfg_file_path.strpath))

  inifile[u"Tabulation"][u'target'] = u"lammps_eam_alloy"

  modified = io.StringIO()
  inifile.write(modified)
  modified.seek(0)

  inifile = configparser.ConfigParser()
  inifile.read_file(modified)
  assert inifile[u"Tabulation"][u'target'] == u"lammps_eam_alloy"

  modified.seek(0)
  config_parser = ConfigParser(modified)
  assert config_parser.tabulation.target == u"setfl"

@needsLAMMPS
def test_lammps_setfl_crg_tabulate_ThO2(lammps_run_fluorite_fixture):
  tmpdir = lammps_run_fluorite_fixture

  cfgobj = Configuration()
  config_file = io.open(_get_lammps_resource_dir().join("CRG_U_Th.aspot").strpath, encoding = "utf8")
  tabulation = cfgobj.read(config_file)

  with tmpdir.join("table.eam.alloy").open("w") as outfile:
    tabulation.write(outfile)

  with lammps_run_fluorite_fixture.join("potentials.lmpinc").open('w') as potfile:
    potfile.write(u"""variable O equal 1
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
  config_file =io.open( _get_lammps_resource_dir().join("CRG_U_Th.aspot").strpath, encoding = "utf8")
  tabulation = cfgobj.read(config_file)

  with tmpdir.join("table.eam.alloy").open("w") as outfile:
    tabulation.write(outfile)

  with lammps_run_fluorite_fixture.join("potentials.lmpinc").open('w') as potfile:
    potfile.write(u"""variable O equal 1
set type 1 charge -1.1104
set type 2 charge 2.2208

kspace_style pppm 1.0e-6

pair_style hybrid/overlay coul/long 10.0 eam/alloy
pair_coeff * * coul/long
pair_coeff * * eam/alloy table.eam.alloy O U
""")

  runLAMMPS(cwd = tmpdir.strpath)
  energy = extractLAMMPSEnergy(cwd = tmpdir.strpath)

  expect = -163.072240194504
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
  config_file = io.open(rd.join("CRG_Ce.aspot").strpath, encoding = "utf8")
  tabulation = cfgobj.read(config_file)

  with tmpdir.join("TABEAM").open("w") as outfile:
    tabulation.write(outfile)

  runDLPoly(cwd = tmpdir.strpath)
  actual = extractDLPOLYEnergy(cwd = tmpdir.strpath)

  expect = -532.6778
  assert pytest.approx(expect) == actual

@needsLAMMPS
def test_lammps_EAM_FS_tabulate_AlFe(lammps_run_fixture):
  tmpdir = lammps_run_fixture

  cfgobj = Configuration()
  config_file = io.open(_get_lammps_resource_dir().join("AlFe_setfl_fs.aspot").strpath, encoding = "utf8")
  tabulation = cfgobj.read(config_file)

  _get_lammps_resource_dir().join("random_Al_Fe.lmpstruct").copy(tmpdir.join("structure.lmpstruct"))
  _get_lammps_resource_dir().join("AlFe_mm.eam.fs").copy(tmpdir.join("table.eam.fs"))

  with tmpdir.join("potentials.lmpinc").open("w") as potfile:
    potfile.write(u"pair_style eam/fs\n")
    potfile.write(u"pair_coeff * * table.eam.fs Al Fe\n")

  runLAMMPS(cwd = tmpdir.strpath)
  expect = extractLAMMPSEnergy(cwd = tmpdir.strpath)

  with tmpdir.join("table.eam.fs").open("w") as outfile:
    tabulation.write(outfile)
  
  runLAMMPS(cwd = tmpdir.strpath)
  actual = extractLAMMPSEnergy(cwd = tmpdir.strpath)

  assert pytest.approx(expect) == actual

@needsDLPOLY
def test_dlpoly_EAM_FS_tabulate_AlFe(tmpdir):
  cfg_file_path = _get_lammps_resource_dir().join("AlFe_setfl_fs.aspot")

  from atsim.potentials.config._config_parser import _RawConfigParser
  inifile = _RawConfigParser()
  inifile.read(u"{}".format(cfg_file_path.strpath))

  inifile[u"Tabulation"][u'target'] = u"DL_POLY_EAM_fs"

  modified = io.StringIO()
  inifile.write(modified)
  modified.seek(0)

  cfgobj = Configuration()
  tabulation = cfgobj.read(modified)

  with tmpdir.join("TABEAM").open("w") as outfile:
    tabulation.write(outfile)

  _get_dlpoly_resource_dir().join("CONTROL_random_Al_Fe").copy(tmpdir.join("CONTROL"))
  _get_dlpoly_resource_dir().join("CONFIG_random_Al_Fe").copy(tmpdir.join("CONFIG"))
  _get_dlpoly_resource_dir().join("FIELD_random_Al_Fe").copy(tmpdir.join("FIELD"))

  runDLPoly(cwd = tmpdir.strpath)
  actual = extractDLPOLYEnergy(cwd = tmpdir.strpath)

  expect = -31.769632

  assert pytest.approx(expect) == actual

def test_custom_species_data():
  cfg_file_path = _get_lammps_resource_dir().join("CRG_U_Th.aspot")

  from atsim.potentials.config._config_parser import _RawConfigParser
  inifile = _RawConfigParser()
  inifile.read(u"{}".format(cfg_file_path.strpath))

  inifile.add_section(u'Species')
  inifile[u"Species"][u"U.atomic_mass"] = u"235"
  inifile[u"Species"][u"U.lattice_constant"] = u"5.678"
  inifile[u"Species"][u"Th.lattice_type"] = u"bcc"

  modified = io.StringIO()
  inifile.write(modified)
  modified.seek(0)

  cfgobj = Configuration()
  tabulation = cfgobj.read(modified)

  plist = tabulation.eam_potentials
  upot = [p for p in plist if p.species == u"U"][0]

  assert pytest.approx(235.0) == upot.mass
  assert pytest.approx(5.678) == upot.latticeConstant

  cepot = [p for p in plist if p.species == u"Th"][0]
  assert cepot.latticeType == u'bcc'
