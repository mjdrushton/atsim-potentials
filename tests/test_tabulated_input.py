import pytest

from io import StringIO

from deepdiff import DeepDiff


example_input_1 = u"""

[Table-Form:cubic_spline]
interpolation : cubic_spline
x : 0 0.01370 0.02740 0.05481 0.06303
y : 0.0 -2.9239 -4.2953 -2.8523 0.0

"""


example_input_2 = u"""

[Table-Form:cubic_spline]
interpolation : cubic_spline
xy : 0.0      0.0
     0.01370 -2.9239
     0.02740 -4.2953
     0.05481 -2.8523
     0.06303 0.0

"""

example_input_3 = u"""

[Table-Form:cubic_spline]
x : 0 0.01370 0.02740 0.05481 0.06303
y : 0.0 -2.9239 -4.2953 -2.8523 0.0

"""

#   * Two table forms with the same name 
bad_input_1 = u"""

[Table-Form:same_name]
x : 1 2 3
y : 1 2 3

[Table-Form:same_name]
x : 1 2 3
y : 1 2 3

"""

#   * Two table forms with the same name (whitespace)
bad_input_2 = u"""

[Table-Form:same_name ]
x : 1 2 3
y : 1 2 3

[Table-Form: same_name]
x : 1 2 3
y : 1 2 3

"""

#   * xy, x and y and combinationas thereof specified at the same time
bad_input_3 = u"""

[Table-Form:bad-input-3]
x : 1 2 3
y : 1 2 3
xy : 0.0      0.0
     0.01370 -2.9239
     0.02740 -4.2953
     0.05481 -2.8523
     0.06303 0.0
"""

bad_input_4 = u"""

[Table-Form:bad-input-4]
x : 1 2 3
xy : 0.0      0.0
     0.01370 -2.9239
     0.02740 -4.2953
     0.05481 -2.8523
     0.06303 0.0
"""

bad_input_5 = u"""

[Table-Form:bad-input-5]
y : 1 2 3
xy : 0.0      0.0
     0.01370 -2.9239
     0.02740 -4.2953
     0.05481 -2.8523
     0.06303 0.0
"""

#   * Different number of points in x and y lists
bad_input_6 = u"""

[Table-Form:bad-input-6]
x : 1 2 
y : 1 2 3

"""

bad_input_7 = u"""

[Table-Form:bad-input-7]
x : 1 2 3
y : 1 2 

"""

#   * Length of xy is not even (isn't pairs)
bad_input_8 = u"""

[Table-Form:bad-input-8]
xy : 0.0      0.0
     0.01370 -2.9239
     0.02740 -4.2953
     0.05481 -2.8523
     0.06303 
"""

bad_input_9 = u"""

[Table-Form:bad-input-9]
x : 1 b 3
y : 1 2 c

"""

# insufficient_length = u"""

# [Table-Form:insufficient-length]
# x : 1 2
# y : 1 2

# """



# example_input_3 = u"""

# [Table-Form:cubic_spline]
# interpolation : cubic_spline
# filename : cubic_spline.csv

# """

# from atsim.potentials.config._table_form_registry import Table_Form_Registry
from atsim.potentials.config import ConfigParser
from atsim.potentials.config._common import  TableFormTuple
from atsim.potentials.config._common import  ConfigParserException, ConfigParserDuplicateEntryException
from atsim.potentials.config._config_parser import _TableFormSection

@pytest.mark.parametrize("cfgstring", [example_input_1, example_input_2, example_input_3])
def test_config_parser(cfgstring):
  cfgio = StringIO(cfgstring)
  cfgparser =  ConfigParser(cfgio)

  expect = [
    TableFormTuple(
      name = u'cubic_spline',
      interpolation = u'cubic_spline',
      x = [0.0, 0.01370, 0.02740, 0.05481, 0.06303],
      y = [0.0, -2.9239, -4.2953, -2.8523, 0.0])
  ]

  actual = cfgparser.table_form
  assert DeepDiff(expect, actual) == {}

def test_orphan_section():
  # Check that Table-Form sections do not appear in config parser's orphan sections.
  cfgio = StringIO(example_input_1)
  cfgparser =  ConfigParser(cfgio)

  for section in cfgparser.orphan_sections:
    if section.startswith(_TableFormSection._section_name_prefix):
      pytest.fail("Section '{}' listed as an orphan section".format(section))

@pytest.mark.parametrize("cfgstring", [bad_input_1, bad_input_2])
def test_duplicate_entry_exception(cfgstring):
  cfgio = StringIO(cfgstring)
  with pytest.raises(ConfigParserDuplicateEntryException):
    ConfigParser(cfgio)

@pytest.mark.parametrize("cfgstring", [bad_input_1, bad_input_2])
def test_bad_input_to_cfgparser_duplicate_sections(cfgstring):
  cfgio = StringIO(cfgstring)
  with pytest.raises(ConfigParserException):
    ConfigParser(cfgio)

@pytest.mark.parametrize("cfgstring", 
  [bad_input_3, bad_input_4, bad_input_5, bad_input_6, bad_input_7, bad_input_8, bad_input_9])
def test_bad_input_to_cfgparser(cfgstring):
  cfgio = StringIO(cfgstring)
  cfgparser =  ConfigParser(cfgio)

  with pytest.raises(ConfigParserException):
    cfgparser.table_form

def test_tableform_builder():
  from atsim.potentials.config._table_form_builder import Table_Form_Builder
  from atsim.potentials.tableforms import Cubic_Spline_Table_Form
  tfb = Table_Form_Builder()
  assert tfb._config_name_to_class('cubic_spline') == Cubic_Spline_Table_Form

  # Convert from TableFormTuple to and Existing_Potential_Form object
  defn_tuple = TableFormTuple('my_table', 'cubic_spline', [0,1,2,3], [0,1,2,3])
  defn = tfb.create_potential_form(defn_tuple)

  assert defn.signature.parameter_names == ["r"]
  assert pytest.approx(1) == defn()(1)

def test_tableforms_in_potentialformregistry():
  from atsim.potentials.config._potential_form_registry import Potential_Form_Registry

  cfg_string = u"""

[Table-Form:tabulated_1]
interpolation : cubic_spline
x : 0.0 0.00728 0.01455 0.02910 0.03347
y : 0.0 -3.2170 -4.6278 -2.7699 0.0

[Table-Form:tabulated_2]
interpolation : cubic_spline
x : 0 1 2 3
y : 0 1 2 3

"""

  sio = StringIO(cfg_string)
  cfg = ConfigParser(sio)
  pfr = Potential_Form_Registry(cfg, register_standard= False)
  assert pfr.registered == ["tabulated_1", "tabulated_2"]
  tabulated_2 = pfr["tabulated_2"]
  assert pytest.approx(tabulated_2()(1.0)) == 1.0

  # Now do a test to make sure other potential forms can use Table-Forms
  cfg_string = u"""

[Potential-Form]
test(r) = tabulated(r) + 5

[Table-Form:tabulated]
interpolation : cubic_spline
x : 0 1 2 3
y : 0 1 2 3

"""

  sio = StringIO(cfg_string)
  cfg = ConfigParser(sio)
  pfr = Potential_Form_Registry(cfg, register_standard= False)

  assert pfr.registered == ["tabulated", "test"]
  assert pytest.approx(pfr["tabulated"]()(1.0)) == 1.0
  assert pytest.approx(pfr["test"]()(1.0)) == 6.0

def test_cubic_spline_table_form():
  # Au data from Foiles1985
  rho = [0.0, 0.00728, 0.01455, 0.02910, 0.03347]
  F_rho = [0.0, -3.2170, -4.6278, -2.7699, 0.0]

  from atsim.potentials.tableforms import Cubic_Spline_Table_Form

  cubic_spline = Cubic_Spline_Table_Form(rho, F_rho)

  # Make sure values are the same at data-points
  for in_rho, expect_F in zip(rho, F_rho):
    assert pytest.approx(expect_F) == cubic_spline(in_rho)

  # Now check that derivatives are correct (deriv and deriv2)
  assert hasattr(cubic_spline, "deriv")

  # ... and do some tests for interpolated points
  line = Cubic_Spline_Table_Form([0,1,2,3], [2,4,6,8])
  assert pytest.approx(4) == line(1)
  assert pytest.approx(3) == line(0.5)
  assert pytest.approx(2) == line.deriv(2.5)
  assert pytest.approx(0) == line.deriv2(2.5)
