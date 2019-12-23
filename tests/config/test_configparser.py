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
from atsim.potentials.config._common import EAMFSDensitySpeciesTuple, PotentialFormTuple, PotentialFormSignatureTuple
from atsim.potentials.config._common import ConfigurationException, ConfigParserDuplicateEntryException

from ._common import _get_lammps_resource_dir, _get_dlpoly_resource_dir


def test_potential_forms():
  """Testing reading of potential forms from [Potential-Form]"""

  cfg_string = u"""[Potential-Form]
buck_morse(r_ij, A,rho,C,D,gamma,r0) : buck(r_ij, A,rho,C) + morse(r_ij, gamma,r0,D)
density(r_ij, n) : (n/r_ij^8) * (1/2)*(1+erf(20*(r_ij-1.5)))
"""

  expect = [
    PotentialFormTuple(PotentialFormSignatureTuple(u"buck_morse", [u"r_ij", u"A", u"rho", u"C", u"D", u"gamma", u"r0"], False), u"buck(r_ij, A,rho,C) + morse(r_ij, gamma,r0,D)"),
    PotentialFormTuple(PotentialFormSignatureTuple(u"density", [u"r_ij", u"n"], False), u"(n/r_ij^8) * (1/2)*(1+erf(20*(r_ij-1.5)))")]

  parsed = ConfigParser(io.StringIO(cfg_string))
  actual = parsed.potential_form
  assert DeepDiff(expect,actual) == {}

def test_parse_potential_form_signature():
  pfstr = "buck_morse(r_ij, A,rho,C,D,gamma,r0)"

  cfg = ConfigParser(io.StringIO())
  expect = PotentialFormSignatureTuple(u"buck_morse", [u"r_ij", u"A", u"rho", u"C", u"D", u"gamma", u"r0"], False)
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
  assert type(parsed.tabulation.nr) is int

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
  parsed = ConfigParser(io.StringIO(u"""
[Species]
Gd.atomic_mass = 1.234
Gd.charge = 4.567

NewSpecies.atomic_number = 600

Al.lattice_type = fcc
Al.lattice_constant = 7.8910
  
  """))

  expect =  { u'Gd' : {u'atomic_mass' : 1.234, u'charge' : 4.567 },
              u'NewSpecies' : {u'atomic_number' : 600},
              u'Al' : {u'lattice_type' : u'fcc', u'lattice_constant' : 7.8910}
  }

  actual = parsed.species

  assert DeepDiff(expect, actual) == {}

  # Test that an exception is raised if invalid syntax used

  parsed = ConfigParser(io.StringIO(u"""
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
  with io.open(_get_lammps_resource_dir().join("AlFe_setfl_fs.aspot").strpath, encoding="utf8") as infile:
    sio = io.StringIO(infile.read())

  sio.seek(0,os.SEEK_END)
  sio.write(u"\n\n[Charges]\nU : 2.22\nO: -1.11\n")
  sio.write(u"[Masses]\nU : 1.234\nO: 4.56\n")
  sio.seek(0)

  expect = [u'tabulation', u'potential_form', u'eam_embed', u'eam_density_fs', u'pair']
  expect.sort()

  cp = ConfigParser(sio)
  actual = cp.parsed_sections
  actual.sort()
  assert expect == actual

  actual = cp.orphan_sections
  expect = [u"Charges", u"Masses"]

  actual.sort()
  expect.sort()

  assert expect == actual


def test_overrides():
  # Test changing a values
  with io.open(_get_lammps_resource_dir().join("zbl_spline.aspot").strpath, encoding = "utf8") as infile:
    cp = ConfigParser(infile, [
      ConfigParserOverrideTuple(u"Tabulation", u"target", u"DLPOLY"),
      ConfigParserOverrideTuple(u"Potential-Form", u"bks(r,qi,qj,A,rho,C)", u"as.buck(r, A, rho, C)")
    ])

  assert "DLPOLY" == cp.tabulation.target
  assert 5000 == cp.tabulation.nr
  assert 10.0 == cp.tabulation.cutoff

  assert 1 == len(cp.potential_form)
  k,v = cp.potential_form[0]

  assert k.label == "bks"
  assert k.parameter_names == [u"r",u"qi",u"qj",u"A",u"rho",u"C"]
  assert v == u"as.buck(r, A, rho, C)"
  
  # Make sure an error is thrown when attempting to change a value that doesn't exist
  with io.open(_get_lammps_resource_dir().join("zbl_spline.aspot").strpath, encoding = "utf8") as infile:
    with pytest.raises(ConfigOverrideException):
      cp = ConfigParser(infile, overrides = [
        ConfigParserOverrideTuple(u"Tabulation", u"targe", u"DLPOLY"),
      ])

  # Test adding extra value to a section
  with io.open(_get_lammps_resource_dir().join("zbl_spline.aspot").strpath, encoding = "utf8") as infile:
      cp = ConfigParser(infile, additional = [
        ConfigParserOverrideTuple(u"Tabulation", u"blah", u"blah"),
      ])

  expect = [u'target', u'cutoff', u'nr', u'blah']
  actual = list(cp.raw_config_parser[u'Tabulation'])
  assert expect == actual
  assert 'blah' == cp.raw_config_parser[u'Tabulation'][u'blah']

  # Make sure an error is thrown if trying to add when item already exists
  with io.open(_get_lammps_resource_dir().join("zbl_spline.aspot").strpath, encoding = "utf8") as infile:
    with pytest.raises(ConfigOverrideException):
      cp = ConfigParser(infile, additional = [
        ConfigParserOverrideTuple(u"Tabulation", u"nr", u"blah"),
      ])

  # test removing a value from a section
  with io.open(_get_lammps_resource_dir().join("zbl_spline.aspot").strpath, encoding = "utf8") as infile:
      cp = ConfigParser(infile, overrides = [
        ConfigParserOverrideTuple(u"Tabulation", u"nr", None),
      ])
      
  with io.open(_get_lammps_resource_dir().join("zbl_spline.aspot").strpath, encoding = "utf8") as infile:
    with pytest.raises(ConfigOverrideException):
      cp = ConfigParser(infile, overrides = [
        ConfigParserOverrideTuple(u"Tabulation", u"no-exist", None),
      ])

  expect = [u'target', u'cutoff', ]
  actual = list(cp.raw_config_parser[u'Tabulation'])
  assert expect == actual

  # test removing the last item from a section
  with io.open(_get_lammps_resource_dir().join("zbl_spline.aspot").strpath, encoding = "utf8") as infile:
      cp = ConfigParser(infile, overrides = [
        ConfigParserOverrideTuple(u"Tabulation", u"nr", None),
        ConfigParserOverrideTuple(u"Tabulation", u"target", None),
        ConfigParserOverrideTuple(u"Tabulation", u"cutoff", None),
      ])

  assert not cp.raw_config_parser.has_section(u"Tabulation")

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

def test_duplicate_entries():
  assert issubclass(ConfigParserDuplicateEntryException, ConfigParserException)
  assert issubclass(ConfigParserDuplicateEntryException, ConfigurationException)

  aspot = u"""

[Pair]
U-U : as.buck 1000.0 0.1 32.0
U-U : as.buck 1000.0 0.2 32.0"""

  with pytest.raises(ConfigParserDuplicateEntryException):
    cp = ConfigParser(io.StringIO(aspot))


  aspot = u"""

[Pair]
O-U : as.buck 1000.0 0.1 32.0
U-O : as.buck 1000.0 0.2 32.0"""

  with pytest.raises(ConfigParserDuplicateEntryException):
    cp = ConfigParser(io.StringIO(aspot))

  aspot = u"""

[Pair]
O-U : as.buck 1000.0 0.1 32.0
U-U : as.buck 2000.0 0.2 16.0"""

  with pytest.raises(ConfigParserDuplicateEntryException):
    cp = ConfigParser(io.StringIO(aspot), additional=[ConfigParserOverrideTuple(u"Pair", u"O-U", u"as.buck 1000.0 0.1 32.0")])

  aspot = u"""

[Pair]
O-U : as.buck 1000.0 0.1 32.0
U-U : as.buck 2000.0 0.2 16.0"""

  with pytest.raises(ConfigParserDuplicateEntryException):
    cp = ConfigParser(io.StringIO(aspot), additional=[ConfigParserOverrideTuple(u"Pair", u"U-O", u"as.buck 1000.0 0.1 32.0")])


def test_indented_input():
  single_line = u"""[Pair]
O-O : sum( as.constant 1.0, as.constant 2.0 )
  """

  expect = ConfigParser(io.StringIO(single_line)).pair

  cfg_string = u"""[Pair]
O-O : sum(
          as.constant
                      1.0,
          as.constant 2.0
      )
  """

  cp = ConfigParser(io.StringIO(cfg_string))
  actual = cp.pair

  assert DeepDiff(expect, actual) == {}

def test_indented_input_with_comment():
  single_line = u"""[Pair]
O-O : sum( as.constant 1.0, as.constant 2.0 )
  """

  expect = ConfigParser(io.StringIO(single_line)).pair
  cfg_string = u"""[Pair]
O-O : sum(
          # Interleaved comment
          as.constant
                      1.0,
          as.constant 2.0
      )
  """

  cp = ConfigParser(io.StringIO(cfg_string))
  actual = cp.pair

  assert DeepDiff(expect, actual) == {}