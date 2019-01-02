from atsim.potentials.config import ConfigParser
from atsim.potentials.config._common import PotentialFormInstanceTuple, \
  MultiRangeDefinitionTuple, \
  PotentialModifierTuple, \
  PairPotentialTuple, \
  EAMEmbedTuple, \
  EAMFSDensitySpeciesTuple, \
  EAMDensityTuple, \
  SpeciesTuple

PMT = PotentialModifierTuple
PPT = PairPotentialTuple
ST = SpeciesTuple  
MRD = MultiRangeDefinitionTuple
PFI = PotentialFormInstanceTuple


from atsim.potentials.config._common import ConfigParserException

import pytest
from deepdiff import DeepDiff

import io

def test_pair_potential_params():
  """Test reading of potential parameters from [Pair] section of config file."""
  cfg_string = u"""[Pair]
A-B : buck 1000.0 0.1 1.0
B-C : buck 2000.0 0.2 2.0
C-D : morse 3000.0 0.3 3.0
"""
  parsed = ConfigParser(io.StringIO(cfg_string))

  expect = [
    PPT(ST(u"A",u"B"), PFI(u"buck", [1000.0, 0.1, 1.0], MRD(u">", 0.0), None)),
    PPT(ST(u"B",u"C"), PFI(u"buck", [2000.0, 0.2, 2.0], MRD(u">", 0.0), None)),
    PPT(ST(u"C",u"D"), PFI(u"morse", [3000.0, 0.3, 3.0], MRD(u">", 0.0), None))]

  actual = parsed.pair
  assert DeepDiff(expect, actual) == {}

def test_eam_embed():
  """Test reading of embedding parameters from [EAM-Embed] section"""
  cfg_string = u"""[EAM-Embed]
A : embed1 1000.0 0.1 1.0
B : embed2 2000.0 0.2 2.0
C = embed3 3000.0 0.3 3.0
"""
  parsed = ConfigParser(io.StringIO(cfg_string))

  expect = [
    EAMEmbedTuple(u"A", PFI(u"embed1", [1000.0, 0.1, 1.0], MRD(u">",0.0), None)),
    EAMEmbedTuple(u"B", PFI(u"embed2", [2000.0, 0.2, 2.0], MRD(u">",0.0), None)),
    EAMEmbedTuple(u"C", PFI(u"embed3", [3000.0, 0.3, 3.0], MRD(u">",0.0), None))]

  actual = parsed.eam_embed
  assert DeepDiff(expect, actual) == {}

def test_eam_density():
  """Test reading of density parameters from [EAM-Density] section"""
  cfg_string = u"""[EAM-Density]
A : density1 1000.0 0.1 1.0
B : density2 2000.0 0.2 2.0
"""
  parsed = ConfigParser(io.StringIO(cfg_string))

  expect = [
    EAMDensityTuple(u"A", PFI(u"density1", [1000.0, 0.1, 1.0], MRD(u">",0.0), None)),
    EAMDensityTuple(u"B", PFI(u"density2", [2000.0, 0.2, 2.0], MRD(u">",0.0), None))]

  actual = parsed.eam_density
  assert DeepDiff(expect, actual) == {}

def test_eam_fs_density():
  """Test reading of density parameter from [EAM-Density] section of config assuming Finnis-Sinclair potential"""

  cfg_string = u"""[EAM-Density]
A->A : density1 1000.0 0.1 1.0
B->C : density2 2000.0 0.2 2.0
"""
  parsed = ConfigParser(io.StringIO(cfg_string))

  expect = [
    EAMDensityTuple(EAMFSDensitySpeciesTuple(u"A", u"A"), PotentialFormInstanceTuple(u"density1", [1000.0, 0.1, 1.0], MultiRangeDefinitionTuple(u">",0.0), None)),
    EAMDensityTuple(EAMFSDensitySpeciesTuple(u"B", u"C"), PotentialFormInstanceTuple(u"density2", [2000.0, 0.2, 2.0], MultiRangeDefinitionTuple(u">",0.0), None))]

  actual = parsed.eam_density_fs
  assert DeepDiff(expect, actual) == {}

def test_multi_range_potential_form():
  """Tests definition of multiple ranges for potential-form definitions"""

  k = u"A"
  v = u"potential 1.0 2.0 3.0"

  parser = ConfigParser(io.StringIO())
  actual = parser._parse_multi_range(k, v)
  assert actual.species == k
  assert actual.potential_form_instance.potential_form == u"potential"
  assert actual.potential_form_instance.parameters==[1.0, 2.0, 3.0]
  assert actual.potential_form_instance.next is None
  assert actual.potential_form_instance.start == (u'>', 0.0)

  k = u"A"
  v = u">=0 potential  1.0 2.0 3.0"

  actual = parser._parse_multi_range(k, v)
  assert actual.species == k
  assert actual.potential_form_instance.potential_form == u"potential"
  assert actual.potential_form_instance.parameters==[1.0, 2.0, 3.0]
  assert actual.potential_form_instance.next is None
  assert actual.potential_form_instance.start == (u'>=', 0.0)

  k = u"A"
  v = u">=0 potential >10 potentialb 1.0 2.0 3.0"

  actual = parser._parse_multi_range(k, v)
  assert actual.species == k
  assert actual.potential_form_instance.potential_form == u"potential"
  assert actual.potential_form_instance.parameters==[]
  assert actual.potential_form_instance.start == (u'>=', 0.0)
  assert actual.potential_form_instance.next.potential_form == u"potentialb"
  assert actual.potential_form_instance.next.start == (u">", 10.0)
  assert actual.potential_form_instance.next.parameters == [1.0,2.0,3.0]

  v = u">= 0.0 potential 1.0 2.0 3.0"
  actual = parser._parse_multi_range(k, v)
  assert actual.potential_form_instance.start == (u'>=', 0.0)

  k = u"A"
  v = u"potential 1.0 2.0 3.0 >1e1 potentialb 5.0 6.0 7.0 >=2.0E1 zero"

  actual = parser._parse_multi_range(k, v)
  assert actual.species == k
  assert actual.potential_form_instance.potential_form == u"potential"
  assert actual.potential_form_instance.parameters==[1.0, 2.0, 3.0]
  assert not actual.potential_form_instance.next is None
  assert actual.potential_form_instance.start == (u'>', 0.0)

  actual = actual.potential_form_instance.next
  assert actual.potential_form == u"potentialb"
  assert actual.parameters == [5.0, 6.0, 7.0]
  assert not actual.next is None
  assert actual.start == (u'>', 10.0)

  actual = actual.next
  assert actual.potential_form == u"zero"
  assert actual.parameters == []
  assert actual.next is None
  assert actual.start == (u'>=', 20.0)

  k = u"A"
  v = u">0.01 potential 1.0"
  actual = parser._parse_multi_range(k, v)
  assert actual.potential_form_instance.start == (u'>', 0.01)
  v = u">1e-2 potential 1.0"
  actual = parser._parse_multi_range(k, v)
  assert actual.potential_form_instance.start == (u'>', 0.01)
  v = u">1.0E-2 potential 1.0"
  actual = parser._parse_multi_range(k, v)
  assert actual.potential_form_instance.start == (u'>', 0.01)
  # v = ">.1E-1 potential 1.0"
  # actual = parser._parse_multi_range(k, v)
  # assert actual.potential_form_instance.start == ('>', 0.01)

def test_sum_modifier():
  """Test instantiation of potential forms that use the sum() modifier and multiple ranges"""
  k = u"A"
  v = u"sum( as.constant 1.0 >=1.0 as.constant 2.0, >1.5 as.constant 3.0)"
  
  parser = ConfigParser(io.StringIO())

  potential_forms = [
    PotentialFormInstanceTuple(potential_form = u'as.constant', parameters = [1.0], start = MultiRangeDefinitionTuple(range_type = u">", start = 0.0),
      next = PotentialFormInstanceTuple(potential_form = u'as.constant', parameters = [2.0], start = MultiRangeDefinitionTuple(range_type = u">=", start = 1.0),
              next = None)),
    PotentialFormInstanceTuple(potential_form = u'as.constant', parameters = [3.0], start = MultiRangeDefinitionTuple(range_type = u">", start = 1.5), next = None) ]

  expect = PotentialModifierTuple(modifier = u"sum", potential_forms = potential_forms,
    start = MultiRangeDefinitionTuple(range_type = u">", start = 0.0),
    next = None)

  expect = PairPotentialTuple(species = k, potential_form_instance = expect)
  actual = parser._parse_multi_range(k, v)
  assert DeepDiff(expect, actual) == {}

def test_sum_modifier_as_part_of_multi_range():
  """Test instantiation of potential forms that use the sum() modifier as sub-range in potential definition"""
  k = u"A"
  v = u"as.constant 2.0 >=1.0 sum(as.constant 1.0 >= 2.0 as.constant 0.5, >=1.5 as.constant 10.0) >= 3.0 zero"
  
  parser = ConfigParser(io.StringIO())

  potential_forms = [
    PFI(u'as.constant',[1.0], MRD(u'>', 0.0),
      PFI(u'as.constant', [0.5], MRD(u'>=', 2.0), None)),
    PFI(u'as.constant', [10.0], MRD(u'>=', 1.5), None)]

  expect = PFI(u'as.constant',[2.0], MRD(u'>', 0.0), 
    PMT(u'sum', potential_forms, MRD(u'>=', 1.0),
      PFI(u'zero', [], MRD(u'>=', 3.0), None)))

  expect = PairPotentialTuple(species = k, potential_form_instance = expect)
  actual = parser._parse_multi_range(k, v)
  assert DeepDiff(expect, actual) == {}

def test_modifier_with_nested_modifier():
  """Test instantiation of potential forms that has a modifier that has a nested modifier in its argument list"""
  k = u"A"
  v = u"zero >=1.0 sum(nested(constant 1.0 >2 zero), buck 10.0 0.1 32.0)"
  
  parser = ConfigParser(io.StringIO())

  potential_forms = [
    PMT(u'nested', [PFI(u'constant', [1.0], MRD(u'>',0.0),
                        PFI(u'zero', [], MRD(u'>', 2), None))], MRD(u'>', 0.0), None),
    PFI(u'buck', [10.0, 0.1, 32.0], MRD(u'>', 0.0), None)]

  expect = PFI(u'zero', [], MRD(u'>', 0.0),
               PMT(u'sum', potential_forms, MRD(u'>=', 1.0), None))

  expect = PairPotentialTuple(species = k, potential_form_instance = expect)
  actual = parser._parse_multi_range(k, v)
  assert DeepDiff(expect, actual) == {}

def test_sum_modifier_fs_start():
  k = u"Al->Al"
  v = u">=0 sum(dens4 0.00019850823042883 2.5 >2.5 as.zero)"
  
  parser = ConfigParser(io.StringIO())
  actual = parser._parse_multi_range(k,v)
  assert DeepDiff(MultiRangeDefinitionTuple(u">=",0), actual.potential_form_instance.start) == {}

def test_modifier_parse_exceptions():
  """Check that a ConfigException is raised when parse errors are encountered"""

  parser = ConfigParser(io.StringIO())
  with pytest.raises(ConfigParserException):
    parser._parse_multi_range(u"A", u"potential1 1.0 2.0 3.0 potential2 2.0")
