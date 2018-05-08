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
    (("A","B"), "buck", [1000.0, 0.1, 1.0]),
    (("B","C"), "buck", [2000.0, 0.2, 2.0]),
    (("C","D"), "morse", [3000.0, 0.3, 3.0])]

  actual = parsed.pair
  assert expect == actual

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
