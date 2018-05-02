import io
import pytest

from atsim.potentials._config import ConfigParser

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
