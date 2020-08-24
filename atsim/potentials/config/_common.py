import collections
import sys

try:
  from inspect import signature, Parameter
except ImportError:
  from funcsigs import signature, Parameter

SpeciesTuple = collections.namedtuple("SpeciesTuple", ["species_a", "species_b"])
EAMFSDensitySpeciesTuple = collections.namedtuple("EAMFSDensitySpeciesTuple", ["from_species", "to_species"])

_tuple_args = ["species", 'potential_form_instance']

EAMEmbedTuple = collections.namedtuple("EAMEmbedTuple", _tuple_args)
EAMDensityTuple = collections.namedtuple("EAMDensityTuple", _tuple_args)
PairPotentialTuple = collections.namedtuple("PairPotentialTuple", _tuple_args)

PotentialFormInstanceTuple = collections.namedtuple('PotentialFormInstanceTuple', ['potential_form', 'parameters', 'start', 'next'])
PotentialFormSignatureTuple = collections.namedtuple("PotentialFormSignatureTuple", ["label", "parameter_names", "is_varargs"])
PotentialFormTuple = collections.namedtuple("PotentialFormTuple", ["signature", "expression"])

MultiRangeDefinitionTuple = collections.namedtuple("MultiRangeDefinition", ["range_type", "start"])
PotentialModifierTuple = collections.namedtuple("PotentialModifierTuple", ['modifier', 'potential_forms', 'start', 'next'])

TableFormTuple = collections.namedtuple("TableFormTuple", ['name', 'interpolation', 'x', 'y'])

class ConfigurationException(Exception):
  pass

class ConfigParserException(ConfigurationException):
  pass

class ConfigParserMissingSectionException(ConfigParserException):
  pass

class ConfigParserDuplicateEntryException(ConfigParserException):
  pass

class Potential_Form_Registry_Exception(ConfigurationException):
  pass

class Potential_Form_Exception(ConfigurationException):
  pass

class Table_Form_Exception(ConfigurationException):
  pass

class Modifier_Exception(ConfigurationException):
  pass

class Unknown_Modifier_Exception(Modifier_Exception):
  pass

def _is_vararg_signature(sig):
  for p in sig.parameters.values():
    if not p.kind == Parameter.VAR_POSITIONAL:
      return False
  return True

def make_potential_form_tuple_from_function(name, pyfunc):
  sig = signature(pyfunc)
  
  # Is this a varargs function?
  varargs = _is_vararg_signature(sig)
  args = []

  if not varargs:
    for param in sig.parameters.values():
      args.append(param.name)

  d = PotentialFormTuple(signature = PotentialFormSignatureTuple(name, args, varargs), expression = "")

  return d

