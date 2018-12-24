import pytest

from atsim.potentials import potentialforms
import atsim.potentials

from atsim.potentials.config._common import PotentialModifierTuple, MultiRangeDefinitionTuple, PotentialFormInstanceTuple
from atsim.potentials.config._common import ConfigurationException
from atsim.potentials.config._potential_form_builder import Potential_Form_Builder
from atsim.potentials.config._potential_form_registry import Potential_Form_Registry
from atsim.potentials.config._config_parser import ConfigParser
from atsim.potentials.config._modifier_registry import Modifier_Registry
from atsim.potentials._modifiers import spline

import io

def test_spline_modifier():
  bks_buck = potentialforms.buck(18003.7572, 1.0/4.87318, 133.5381)
  bks_coul = potentialforms.coul(2.4, -1.2)
  bks = atsim.potentials.plus(bks_buck, bks_coul)
  zbl = potentialforms.zbl(14, 8)
  spline = atsim.potentials.SplinePotential( zbl, bks_buck, 0.8, 1.4)
  
  buck_spline = PotentialFormInstanceTuple(potential_form = "as.buck", parameters = [18003.7572, 1.0/4.87318, 133.5381], start = MultiRangeDefinitionTuple('>', 1.4), next = None)
  exp_spline = PotentialFormInstanceTuple(potential_form = "exp_spline", parameters = [], start = MultiRangeDefinitionTuple(">=", 0.8), next = buck_spline)
  zbl_spline = PotentialFormInstanceTuple(potential_form = "as.zbl", parameters = [14, 8], start = MultiRangeDefinitionTuple(">", 0.0), next = exp_spline)

  pmt = PotentialModifierTuple(modifier = 'spline', 
                              potential_forms = [zbl_spline],
                              start = MultiRangeDefinitionTuple('>', 0.0),
                              next = None)

  pfr = Potential_Form_Registry(ConfigParser(io.StringIO()), True)
  mr = Modifier_Registry()
  pfb = Potential_Form_Builder(pfr, mr)
  mod_spline = pfb.create_potential_function(pmt)

  for i in range(1,100):
    r = float(i) / 10.0

    expect = spline(r)
    actual = mod_spline(r)

    assert pytest.approx(expect, abs = 1e-3) == actual

  assert pytest.approx(zbl(0.7)) == mod_spline(0.7)
  assert pytest.approx(bks_buck(2.0)) == mod_spline(2.0)
  assert pytest.approx(spline(1.2)) == mod_spline(1.2)

  assert pytest.approx(zbl.deriv(0.7)) == mod_spline.deriv(0.7)
  assert pytest.approx(bks_buck.deriv(2.0)) == mod_spline.deriv(2.0)
  assert pytest.approx(spline.deriv(1.2)) == mod_spline.deriv(1.2)

  assert pytest.approx(zbl.deriv2(0.7)) == mod_spline.deriv2(0.7)
  assert pytest.approx(bks_buck.deriv2(2.0)) == mod_spline.deriv2(2.0)
  assert pytest.approx(spline.deriv2(1.2)) == mod_spline.deriv2(1.2)
  

def test_bad_ranges():
  # Bad spline name
  buck_spline = PotentialFormInstanceTuple(potential_form = "as.buck", parameters = [18003.7572, 1.0/4.87318, 133.5381], start = MultiRangeDefinitionTuple('>', 1.4), next = None)
  exp_spline = PotentialFormInstanceTuple(potential_form = "bad_spline", parameters = [], start = MultiRangeDefinitionTuple(">=", 0.8), next = buck_spline)
  zbl_spline = PotentialFormInstanceTuple(potential_form = "as.zbl", parameters = [14, 8], start = MultiRangeDefinitionTuple(">", 0.0), next = exp_spline)

  pmt = PotentialModifierTuple(modifier = 'spline', 
                            potential_forms = [zbl_spline],
                            start = MultiRangeDefinitionTuple('>', 0.0),
                            next = None)

  pfr = Potential_Form_Registry(ConfigParser(io.StringIO()), True)
  mr = Modifier_Registry()
  pfb = Potential_Form_Builder(pfr, mr)

  with pytest.raises(ConfigurationException):
    pfb.create_potential_function(pmt)  

  # Too few ranges
  exp_spline = PotentialFormInstanceTuple(potential_form = "exp_spline", parameters = [], start = MultiRangeDefinitionTuple(">=", 0.8), next = None)
  zbl_spline = PotentialFormInstanceTuple(potential_form = "as.zbl", parameters = [14, 8], start = MultiRangeDefinitionTuple(">", 0.0), next = exp_spline)

  pmt = PotentialModifierTuple(modifier = 'spline', 
                            potential_forms = [zbl_spline],
                            start = MultiRangeDefinitionTuple('>', 0.0),
                            next = None)

  with pytest.raises(ConfigurationException):
    pfb.create_potential_function(pmt)  

  # Too many
  buck_spline2 = PotentialFormInstanceTuple(potential_form = "as.buck", parameters = [18003.7572, 1.0/4.87318, 133.5381], start = MultiRangeDefinitionTuple('>', 1.4), next = None)
  buck_spline = PotentialFormInstanceTuple(potential_form = "as.buck", parameters = [18003.7572, 1.0/4.87318, 133.5381], start = MultiRangeDefinitionTuple('>', 1.4), next = buck_spline2)
  exp_spline = PotentialFormInstanceTuple(potential_form = "bad_spline", parameters = [], start = MultiRangeDefinitionTuple(">=", 0.8), next = buck_spline)
  zbl_spline = PotentialFormInstanceTuple(potential_form = "as.zbl", parameters = [14, 8], start = MultiRangeDefinitionTuple(">", 0.0), next = exp_spline)

  pmt = PotentialModifierTuple(modifier = 'spline', 
                            potential_forms = [zbl_spline],
                            start = MultiRangeDefinitionTuple('>', 0.0),
                            next = None)

  pfr = Potential_Form_Registry(ConfigParser(io.StringIO()), True)
  mr = Modifier_Registry()
  pfb = Potential_Form_Builder(pfr, mr)

  with pytest.raises(ConfigurationException):
    pfb.create_potential_function(pmt)  

def test_too_many_potential_forms():
  pfr = Potential_Form_Registry(ConfigParser(io.StringIO()), True)
  mr = Modifier_Registry()
  pfb = Potential_Form_Builder(pfr, mr)
  with pytest.raises(ConfigurationException):
    spline(["one", "two", "three"], pfb)

  with pytest.raises(ConfigurationException):
    spline([], pfb)

def test_buck4_spline():
  pytest.fail()