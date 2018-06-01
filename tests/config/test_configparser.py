import io
import pytest

from atsim.potentials.config import ConfigParser
from atsim.potentials.config import ConfigParserException

def test_pair_potential_params():
  """Test reading of potential parameters from [Pair] section of config file."""
  cfg_string = u"""[Pair]
A-B : buck 1000.0 0.1 1.0
B-C : buck 2000.0 0.2 2.0
C-D : morse 3000.0 0.3 3.0
"""
  parsed = ConfigParser(io.StringIO(cfg_string))

  expect = [
    (("A","B"), ("buck", [1000.0, 0.1, 1.0], (">", 0), None)),
    (("B","C"), ("buck", [2000.0, 0.2, 2.0], (">", 0), None)),
    (("C","D"), ("morse", [3000.0, 0.3, 3.0], (">", 0), None))]

  actual = parsed.pair
  assert expect == actual

def test_eam_embed():
  """Test reading of embedding parameters from [EAM-Embed] section"""
  cfg_string = u"""[EAM-Embed]
A : embed1 1000.0 0.1 1.0
B : embed2 2000.0 0.2 2.0
C = embed3 3000.0 0.3 3.0
"""
  parsed = ConfigParser(io.StringIO(cfg_string))

  expect = [
    ("A", ("embed1", [1000.0, 0.1, 1.0], (">",0), None)),
    ("B", ("embed2", [2000.0, 0.2, 2.0], (">",0), None)),
    ("C", ("embed3", [3000.0, 0.3, 3.0], (">",0), None))]

  actual = parsed.eam_embed
  assert expect == actual

def test_eam_density():
  """Test reading of density parameters from [EAM-Density] section"""
  cfg_string = u"""[EAM-Density]
A : density1 1000.0 0.1 1.0
B : density2 2000.0 0.2 2.0
"""
  parsed = ConfigParser(io.StringIO(cfg_string))

  expect = [
    ("A", ("density1", [1000.0, 0.1, 1.0], (">",0), None)),
    ("B", ("density2", [2000.0, 0.2, 2.0], (">",0), None))]

  actual = parsed.eam_density
  assert expect == actual

def test_multi_range_potential_form():
  """Tests definition of multiple ranges for potential-form definitions"""

  k = "A"
  v = "potential 1.0 2.0 3.0"

  parser = ConfigParser(io.StringIO())
  actual = parser._parse_multi_range(k, v)
  assert actual.species == k
  assert actual.potential_form_instance.potential_form == "potential"
  assert actual.potential_form_instance.parameters==[1.0, 2.0, 3.0]
  assert actual.potential_form_instance.next is None
  assert actual.potential_form_instance.start == ('>', 0.0)

  k = "A"
  v = ">=0 potential  1.0 2.0 3.0"

  actual = parser._parse_multi_range(k, v)
  assert actual.species == k
  assert actual.potential_form_instance.potential_form == "potential"
  assert actual.potential_form_instance.parameters==[1.0, 2.0, 3.0]
  assert actual.potential_form_instance.next is None
  assert actual.potential_form_instance.start == ('>=', 0.0)

  k = "A"
  v = ">=0 potential >10 potentialb 1.0 2.0 3.0"

  actual = parser._parse_multi_range(k, v)
  assert actual.species == k
  assert actual.potential_form_instance.potential_form == "potential"
  assert actual.potential_form_instance.parameters==[]
  assert actual.potential_form_instance.start == ('>=', 0.0)
  assert actual.potential_form_instance.next.potential_form == "potentialb"
  assert actual.potential_form_instance.next.start == (">", 10.0)
  assert actual.potential_form_instance.next.parameters == [1.0,2.0,3.0]

  v = ">= 0.0 potential 1.0 2.0 3.0"
  actual = parser._parse_multi_range(k, v)
  assert actual.potential_form_instance.start == ('>=', 0.0)

  k = "A"
  v = "potential 1.0 2.0 3.0 >1e1 potentialb 5.0 6.0 7.0 >=2.0E1 zero"

  actual = parser._parse_multi_range(k, v)
  assert actual.species == k
  assert actual.potential_form_instance.potential_form == "potential"
  assert actual.potential_form_instance.parameters==[1.0, 2.0, 3.0]
  assert not actual.potential_form_instance.next is None
  assert actual.potential_form_instance.start == ('>', 0.0)

  actual = actual.potential_form_instance.next
  assert actual.potential_form == "potentialb"
  assert actual.parameters == [5.0, 6.0, 7.0]
  assert not actual.next is None
  assert actual.start == ('>', 10.0)

  actual = actual.next
  assert actual.potential_form == "zero"
  assert actual.parameters == []
  assert actual.next is None
  assert actual.start == ('>=', 20.0)

  k = "A"
  v = ">0.01 potential 1.0"
  actual = parser._parse_multi_range(k, v)
  assert actual.potential_form_instance.start == ('>', 0.01)
  v = ">1e-2 potential 1.0"
  actual = parser._parse_multi_range(k, v)
  assert actual.potential_form_instance.start == ('>', 0.01)
  v = ">1.0E-2 potential 1.0"
  actual = parser._parse_multi_range(k, v)
  assert actual.potential_form_instance.start == ('>', 0.01)
  v = ">.1E-1 potential 1.0"
  actual = parser._parse_multi_range(k, v)
  assert actual.potential_form_instance.start == ('>', 0.01)


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
