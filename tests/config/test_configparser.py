import io
import os

import pytest
from deepdiff import DeepDiff

from atsim.potentials.config import ConfigParser, FilteredConfigParser
from atsim.potentials.config import ConfigParserException
from atsim.potentials.config import ConfigParserOverrideTuple
from atsim.potentials.config import ConfigOverrideException
from atsim.potentials.config._config_parser import _RawConfigParser
from atsim.potentials.config._common import SpeciesTuple
from atsim.potentials.config._common import EAMFSDensitySpeciesTuple

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

def test_species():

  # Test empty reference data section
  parsed = ConfigParser(io.StringIO())
  assert parsed.species == {}

  # Test reference data with overrides
  parsed = ConfigParser(io.StringIO("""
[Species]
Gd.atomic_mass = 1.234
Gd.charge = 4.567

NewSpecies.atomic_number = 600

Al.lattice_type = fcc
Al.lattice_constant = 7.8910
  
  """))

  expect =  { 'Gd' : {'atomic_mass' : 1.234, 'charge' : 4.567 },
              'NewSpecies' : {'atomic_number' : 600},
              'Al' : {'lattice_type' : 'fcc', 'lattice_constant' : 7.8910}
  }

  actual = parsed.species

  assert DeepDiff(expect, actual) == {}

  # Test that an exception is raised if invalid syntax used

  parsed = ConfigParser(io.StringIO("""
[Species]
Gd : 1.234
"""))
  with pytest.raises(ConfigParserException):
    parsed.species
    




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


def test_overrides():
  # Test changing a values
  with _get_lammps_resource_dir().join("zbl_spline.aspot").open() as infile:
    cp = ConfigParser(infile, [
      ConfigParserOverrideTuple("Tabulation", "target", "DLPOLY"),
      ConfigParserOverrideTuple("Potential-Form", "bks(r,qi,qj,A,rho,C)", "as.buck(r, A, rho, C)")
    ])

  assert "DLPOLY" == cp.tabulation.target
  assert 5000 == cp.tabulation.nr
  assert 10.0 == cp.tabulation.cutoff

  assert 1 == len(cp.potential_form)
  k,v = cp.potential_form[0]

  assert k.label == "bks"
  assert k.parameter_names == ["r","qi","qj","A","rho","C"]
  assert v == "as.buck(r, A, rho, C)"
  
  # Make sure an error is thrown when attempting to change a value that doesn't exist
  with _get_lammps_resource_dir().join("zbl_spline.aspot").open() as infile:
    with pytest.raises(ConfigOverrideException):
      cp = ConfigParser(infile, overrides = [
        ConfigParserOverrideTuple("Tabulation", "targe", "DLPOLY"),
      ])

  # Test adding extra value to a section
  with _get_lammps_resource_dir().join("zbl_spline.aspot").open() as infile:
      cp = ConfigParser(infile, additional = [
        ConfigParserOverrideTuple("Tabulation", "blah", "blah"),
      ])

  expect = ['target', 'cutoff', 'nr', 'blah']
  actual = list(cp.raw_config_parser['Tabulation'])
  assert expect == actual
  assert 'blah' == cp.raw_config_parser['Tabulation']['blah']

  # Make sure an error is thrown if trying to add when item already exists
  with _get_lammps_resource_dir().join("zbl_spline.aspot").open() as infile:
    with pytest.raises(ConfigOverrideException):
      cp = ConfigParser(infile, additional = [
        ConfigParserOverrideTuple("Tabulation", "nr", "blah"),
      ])


  # test removing a value from a section
  with _get_lammps_resource_dir().join("zbl_spline.aspot").open() as infile:
      cp = ConfigParser(infile, overrides = [
        ConfigParserOverrideTuple("Tabulation", "nr", None),
      ])

  expect = ['target', 'cutoff', ]
  actual = list(cp.raw_config_parser['Tabulation'])
  assert expect == actual

  # test removing the last item from a section
  with _get_lammps_resource_dir().join("zbl_spline.aspot").open() as infile:
      cp = ConfigParser(infile, overrides = [
        ConfigParserOverrideTuple("Tabulation", "nr", None),
        ConfigParserOverrideTuple("Tabulation", "target", None),
        ConfigParserOverrideTuple("Tabulation", "cutoff", None),
      ])

  assert not cp.raw_config_parser.has_section("Tabulation")

def u_th_filtered_test(orig_cp, filtered_cp):
  # First check the unfiltered properties
  unfiltered = [
    "potential_form", 
    "tabulation", 
    "parsed_sections", 
    "orphan_sections", 
    "raw_config_parser"
  ]

  for attr in unfiltered:
    expect = getattr(orig_cp, attr)
    actual = getattr(filtered_cp, attr)
    assert DeepDiff(expect, actual) == {}

  # Now check the attributes that should be filtered
  # ... pair
  ST = SpeciesTuple
  expect = [
    ST("O", "O"),
    ST("U", "U"),
    ST("U", "O") ]
  actual = [p[0] for p in filtered_cp.pair]
  assert DeepDiff(expect, actual) == {}

  # ... eam_embed
  expect = [ "U", "O"]
  actual = [p[0] for p in filtered_cp.eam_embed]  
  assert DeepDiff(expect, actual) == {}

  # ... eam_density
  expect = ["U", "O"]
  actual = [p[0] for p in filtered_cp.eam_density]
  assert DeepDiff(expect, actual) == {}


def test_filtered_config_parser():
  lmpdir = _get_lammps_resource_dir()
  with lmpdir.join("CRG_U_Th.aspot").open() as infile:
    orig_cp = ConfigParser(infile)

  filtered_cp = FilteredConfigParser(orig_cp, exclude = ["Th"])
  u_th_filtered_test(orig_cp, filtered_cp)
  
  filtered_cp = FilteredConfigParser(orig_cp, include = ["U", "O"])
  u_th_filtered_test(orig_cp, filtered_cp)

  with pytest.raises(ValueError):
    FilteredConfigParser(orig_cp, include = ["U", "O"], exclude = ["Th"])


def test_filtered_config_parser_finnis_sinclair():
  lmpdir = _get_lammps_resource_dir()
  with lmpdir.join("AlFe_setfl_fs.aspot").open() as infile:
    orig_cp = ConfigParser(infile)

  filtered_cp = FilteredConfigParser(orig_cp, exclude = ["Fe"])

  unfiltered = [
    "potential_form", 
    "tabulation", 
    "parsed_sections", 
    "orphan_sections", 
    "raw_config_parser"
  ]

  for attr in unfiltered:
    expect = getattr(orig_cp, attr)
    actual = getattr(filtered_cp, attr)
    assert DeepDiff(expect, actual) == {}

  # Now check the attributes that should be filtered
  # ... pair
  ST = SpeciesTuple
  expect = [
    ST("Al", "Al")]
  actual = [p[0] for p in filtered_cp.pair]
  assert DeepDiff(expect, actual) == {}

  # ... eam_embed
  expect = ["Al"]
  actual = [p[0] for p in filtered_cp.eam_embed]  
  assert DeepDiff(expect, actual) == {}

  # ... eam_density
  EDST = EAMFSDensitySpeciesTuple
  expect = [
    EDST("Al", "Al") ]
  actual = [p[0] for p in filtered_cp.eam_density_fs]
  assert DeepDiff(expect, actual) == {}
