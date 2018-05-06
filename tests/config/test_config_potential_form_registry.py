import pytest

import io

from atsim.potentials._config import Potential_Form_Registry, ConfigParser, Potential_Form
from atsim.potentials._config import PotentialFormSignatureTuple, PotentialFormTuple 
from atsim.potentials._config import Potential_Form_Exception, Potential_Form_Registry_Exception

import cexprtk

def test_config_potential_form_registry_one_potential_form():
  cfg_string = u"""[Potential-Form]
buck(r, A, rho, C) : A*exp(-r/rho) - C/r^6
"""
  cfg = ConfigParser(io.StringIO(cfg_string))
  pfr = Potential_Form_Registry(cfg)

  assert ["buck"] == pfr.registered

  buck = pfr["buck"](1000.0, 0.1, 3.0)
  assert pytest.approx(-2.9546) == buck(1.0)


def test_config_potential_form_registry():
  cfg_string = u"""[Potential-Form]
buck(r, A, rho, C) : A*exp(-r/rho) - C/r^6
morse(r, gamma, r_star, D) : D*(exp(-2.0*gamma*(r-r_star)) - 2.0*exp(-gamma*(r-r_star)))
buck_morse(r, A, rho, C, gamma, r_star, D) : buck(r,A,rho,C) + morse(r, gamma, r_star, D)
"""
  cfg = ConfigParser(io.StringIO(cfg_string))
  pfr = Potential_Form_Registry(cfg)

  assert ["buck", "morse", "buck_morse"] == pfr.registered

  buck = pfr["buck"](1000.0, 0.1, 3.0)
  assert pytest.approx(-2.9546) == buck(1.0)

  morse = pfr["morse"](0.719250701502, 1.86874578949, 2.35603812582)
  assert pytest.approx(-0.5811129916666448) == morse(1.0)

  buck_morse = pfr["buck_morse"](1000.0, 0.1, 3.0, 0.719250701502, 1.86874578949, 2.35603812582)
  assert pytest.approx(-0.5811129916666448 + -2.9546) == buck_morse(1.0)

def test_repeat_potential_form_error():
  """Make sure error is raised if two potentials have same label"""
  cfg_string = u"""[Potential-Form]
buck(r, A, rho, C) : A*exp(-r/rho) - C/r^6
buck(r, A, rho, C,D) : A*exp(-r/rho) - C/r^6
"""
  cfg = ConfigParser(io.StringIO(cfg_string))

  with pytest.raises(Potential_Form_Registry_Exception):
    pfr = Potential_Form_Registry(cfg)

def test_potential_instantiation_with_wrong_params_error():
  sig = PotentialFormSignatureTuple('buck', ["r", "A", "rho", "C"])
  pot = PotentialFormTuple(sig, "A*exp(-r/rho) - C/r^6")

  pf = Potential_Form(pot)
  potfunc = pf(1000.0, 0.1, 2.0)
  with pytest.raises(Potential_Form_Exception):
    pf(1000.0, 0.1, 2.0, 10.0)

  with pytest.raises(Potential_Form_Exception):
    pf(2.0, 10.0)



