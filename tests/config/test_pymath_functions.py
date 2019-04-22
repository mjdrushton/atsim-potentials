from atsim.potentials.config import Configuration, ConfigParser
from atsim.potentials.config._potential_form_registry import Potential_Form_Registry

import io
import math

import pytest

base_template = u"""

[Pair]
O-O : test_form

[Potential-Form]
test_form(r_ij) = {}

"""

test_inputs = [
  ("pymath.ceil(4.2)", 5),
  ("pymath.copysign(3, -10)", -3),
  ("pymath.fabs(-3.3)", 3.3),
  ("pymath.factorial(4)", 24),
  ("pymath.floor(-3.1)", -4),
  ("pymath.fmod(2.1, 2)", 0.1),
  # ("pymath.frexp(2.1, 2)", 0.1),
  ("pymath.fsum(1,2,3,4)", 10.0),
  ("pymath.gcd(100,10)", 10.0),
  ("pymath.ldexp(3.1, 4)", 3.1*(2**4)),
  ("pymath.trunc(5.2)", 5.0),
  ("pymath.exp(5.2)", 181.27224187515122),
  ("pymath.log(10)", 2.302585092994046),
  ("pymath.log(10, 10)", 1.0),
  ("pymath.log1p(10)", math.log1p(10.0)),
  ("pymath.log2(10)", 3.321928094887362),
  ("pymath.log2(10)", 3.321928094887362),
  ("pymath.log10(10)", 1.0),
  ("pymath.pow(4,2)", 16.0),
  ("pymath.sqrt(4)", 2.0),
  ("pymath.acos(0.9)", math.acos(0.9)),
  ("pymath.atan(0.9)", math.atan(0.9)),
  ("pymath.atan2(5,2)", math.atan2(5,2)),
  ("pymath.cos(2)", math.cos(2)),
  ("pymath.hypot(3,4)", 5),
  ("pymath.sin(3)", math.sin(3)),
  ("pymath.tan(0.5)", math.tan(0.5)),
  ("pymath.radians(180)", math.radians(180)),
  ("pymath.degrees(3.14)", math.degrees(3.14)),
  ("pymath.acosh(2.9)", math.acosh(2.9)),
  ("pymath.asinh(2.9)", math.asinh(2.9)),
  ("pymath.atanh(0.9)", math.atanh(0.9)),
  ("pymath.cosh(2.9)", math.cosh(2.9)),
  ("pymath.sinh(2.9)", math.sinh(2.9)),
  ("pymath.tanh(0.9)", math.tanh(0.9))
]


@pytest.mark.parametrize("potential_form, expected_output", test_inputs)
def test_pymath(potential_form, expected_output):
  cfg_string = base_template.format(potential_form)
  cfgfile = io.StringIO(cfg_string)

  cfgobj = Configuration()
  tabulation = cfgobj.read(cfgfile)

  pot = tabulation.potentials[0]
  actual = pot.energy(1.0)
  assert pytest.approx(expected_output) == actual

def test_not_registered_as_potentialform():
  cfp = ConfigParser(io.StringIO())
  pfr = Potential_Form_Registry(cfp, register_pymath_functions = True)



