import collections
import sys

SpeciesTuple = collections.namedtuple("SpeciesTuple", ["species_a", "species_b"])
EAMFSDensitySpeciesTuple = collections.namedtuple("EAMFSDensitySpeciesTuple", ["from_species", "to_species"])

_tuple_args = ["species", 'potential_form_instance']
_potential_instantiation_tuples = ['EAMEmbedTuple', 'EAMDensityTuple', 'PairPotentialTuple']

def _populate_module():
  currmodule = sys.modules[__name__]
  for tname in _potential_instantiation_tuples:
    tupletype = collections.namedtuple(tname, _tuple_args)
    setattr(currmodule, tname, tupletype)
_populate_module()

PotentialFormInstanceTuple = collections.namedtuple('PotentialFormInstanceTuple', ['potential_form', 'parameters', 'start', 'next'])
PotentialFormSignatureTuple = collections.namedtuple("PotentialFormSignatureTuple", ["label", "parameter_names", "is_varargs"])
PotentialFormTuple = collections.namedtuple("PotentialFormTuple", ["signature", "expression"])

MultiRangeDefinitionTuple = collections.namedtuple("MultiRangeDefinition", ["range_type", "start"])
PotentialModifierTuple = collections.namedtuple("PotentialModifierTuple", ['modifier', 'potential_forms', 'start', 'next'])

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