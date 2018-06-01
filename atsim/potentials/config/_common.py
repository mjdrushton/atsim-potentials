import collections
import sys

SpeciesTuple = collections.namedtuple("SpeciesTuple", ["species_a", "species_b"])

_tuple_args = ["species", 'potential_form_instance']
_potential_instantiation_tuples = ['EAMEmbedTuple', 'EAMDensityTuple', 'PairPotentialTuple']

def _populate_module():
  currmodule = sys.modules[__name__]
  for tname in _potential_instantiation_tuples:
    tupletype = collections.namedtuple(tname, _tuple_args)
    setattr(currmodule, tname, tupletype)
_populate_module()

PotentialFormInstanceTuple = collections.namedtuple('PotentialFormInstanceTuple', ['potential_form', 'parameters', 'start', 'next'])
PotentialFormSignatureTuple = collections.namedtuple("PotentialFormSignatureTuple", ["label", "parameter_names"])
PotentialFormTuple = collections.namedtuple("PotentialFormTuple", ["signature", "expression"])

MultiRangeDefinitionTuple = collections.namedtuple("MultiRangeDefinition", ["range_type", "start"])

class ConfigurationException(Exception):
  pass

class ConfigParserException(ConfigurationException):
  pass

class ConfigParserMissingSectionException(ConfigParserException):
  pass

class Potential_Form_Registry_Exception(Exception):
  pass

class Potential_Form_Exception(Exception):
  pass