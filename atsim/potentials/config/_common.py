import collections

SpeciesTuple = collections.namedtuple("SpeciesTuple", ["species_a", "species_b"])
PairPotentialTuple = collections.namedtuple("PairPotentialTuple", ["species", "potential_form", "parameters"])
PotentialFormSignatureTuple = collections.namedtuple("PotentialFormSignatureTuple", ["label", "parameter_names"])
PotentialFormTuple = collections.namedtuple("PotentialFormTuple", ["signature", "expression"])

class ConfigParserException(Exception):
  pass

class ConfigParserException(Exception):
  pass

class ConfigParserMissingSectionException(ConfigParserException):
  pass

class Potential_Form_Registry_Exception(Exception):
  pass

class Potential_Form_Exception(Exception):
  pass

PAIR_TABULATION = "Pair"