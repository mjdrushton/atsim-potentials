import pytest

from atsim.potentials.config._potential_form_builder import Potential_Form_Builder
from atsim.potentials.config._potential_form_registry import Potential_Form_Registry
from atsim.potentials.config import ConfigParser
from atsim.potentials.config._common import PotentialFormInstanceTuple, SpeciesTuple
from atsim.potentials.config._common import MultiRangeDefinitionTuple
from atsim.potentials import Multi_Range_Potential_Form

import atsim.potentials.potentialfunctions as pforms


import io

def test_potential_form_builder():
  # Populate potential form registry
  cp = ConfigParser(io.StringIO())
  pfr = Potential_Form_Registry(cp,register_standard = True)

  # Create the potential form builder
  pfb = Potential_Form_Builder(pfr)

  in_tuple = PotentialFormInstanceTuple("as.buck", [1000.0, 0.3, 32.0], None, None)
  potential_func = pfb.create_potential_function(in_tuple)
  assert pytest.approx(pforms.buck(2.0, 1000.0, 0.3, 32.0)) == potential_func(2.0)

def test_multirange_potential_range_search():
  Rdt = Multi_Range_Potential_Form.Range_Defn_Tuple

  tuples = [
  Rdt('>=', 2.0, "one"),
  Rdt('>', float("-inf"), "two"),
  Rdt('>=', -5.0, "three"),
  
  Rdt('>', 3.0, "four"),
  Rdt('>=', 3.0, "five")]

  mrpf = Multi_Range_Potential_Form(*list(tuples))

  assert mrpf._range_search(2.0).potential_form == "one"
  assert mrpf._range_search(-5.0).potential_form == "three"
  assert mrpf._range_search(-100.0).potential_form == "two"
  assert mrpf._range_search(0.0).potential_form == "three"
  assert mrpf._range_search(3.0).potential_form == "five"
  assert mrpf._range_search(3.01).potential_form == "four"

  # Check when no range found for value
  mrpf = Multi_Range_Potential_Form(Rdt(">", 0.0, "one"))
  assert mrpf._range_search(0.0) == None
  assert mrpf._range_search(0.1).potential_form == "one"

  # Check when no tuples are specified
  mrpf = Multi_Range_Potential_Form(default_value = 1.0)
  assert mrpf._range_search(10.0) == None

  tuples = [
    Rdt(range_type='>=', start=0, potential_form="one"), 
    Rdt(range_type='>', start=1, potential_form="two"), 
    Rdt(range_type='>', start=5, potential_form="three")]
  mrpf = Multi_Range_Potential_Form(*tuples)
  actual = mrpf._range_search(1.0)
  assert actual.potential_form == "one"


def test_multirange_potential_set_range_tuples():
  # Check that range tuples appear in ascending order
  Rdt = Multi_Range_Potential_Form.Range_Defn_Tuple
  tuples = [
    Rdt('>=', 2.0, "one"),
    Rdt('>', float("-inf"), "two"),
    Rdt('>=', -5.0, "three"),
    
    Rdt('>', 3.0, "four"),
    Rdt('>=', 3.0, "five")]
  
  expect = ["two", "three", "one", "five", "four"]

  mrpf = Multi_Range_Potential_Form(*list(tuples))
  actual = [rt.potential_form for rt in  mrpf.range_tuples]
  assert expect == actual

  tuples = list(reversed(tuples))
  mrpf.range_tuples = tuples
  actual = [rt.potential_form for rt in  mrpf.range_tuples]
  assert expect == actual
  

def test_multirange_potential_form_builder():
  # Populate potential form registry
  cp = ConfigParser(io.StringIO())
  pfr = Potential_Form_Registry(cp,register_standard = True)

  pfb = Potential_Form_Builder(pfr)

  PFitTup = PotentialFormInstanceTuple
  MRTup = MultiRangeDefinitionTuple

  mrdefn = PFitTup( potential_form = "as.zero", parameters = [], start = MRTup(">=", 0), 
       next = PFitTup( potential_form = "as.buck", parameters = [1000.0, 0.3, 32.0], start = MRTup(">", 1), 
       next = PFitTup( potential_form = "as.constant", parameters = [3], start = MRTup(">", 5), next = None)))

  potential_func = pfb.create_potential_function(mrdefn)
  assert pytest.approx(0) == potential_func(0.0)
  v = potential_func(1.0)
  assert pytest.approx(0) == v
  v = potential_func(2.0)
  assert pytest.approx(pforms.buck(2.0, 1000.0, 0.3, 32.0)) == v
  v = potential_func(5.0)
  assert pytest.approx(pforms.buck(5.0, 1000.0, 0.3, 32.0)) == v
  assert 3.0 == potential_func(5.1)
  assert 3.0 == potential_func(10.0)
