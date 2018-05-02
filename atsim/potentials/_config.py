import collections
import configparser
import re

SpeciesTuple = collections.namedtuple("SpeciesTuple", ["species_a", "species_b"])
PairPotentialTuple = collections.namedtuple("PairPotentialTuple", ["species", "potential_form", "parameters"])
PotentialFormSignatureTuple = collections.namedtuple("PotentialFormSignatureTuple", ["label", "parameter_names"])
PotentialFormTuple = collections.namedtuple("PotentialFormTuple", ["signature", "expression"])

class _RawConfigParser(configparser.RawConfigParser):

  # def __init__(self):
  #   super().__init__()

  def optionxform(self, option):
    return option

class ConfigParserException(Exception):
  pass

class ConfigParser(object):
  """Performs initial stage (tokenizing) of generating a potential model
  suitable for tabulation functions."""

  _signature_re = re.compile(r"^([a-zA-Z]\w*?)\((.*)\)")

  def __init__(self, fp):
    """Construct ConfigParser.

    :param fp: file object containing configuration data."""
    self._config_parser = self._init_config_parser(fp)

  def _init_config_parser(self, fp):
    cp = _RawConfigParser()
    cp.readfp(fp)
    return cp

  def _parse_pair_line(self, k, value):
    species_a, species_b = k.split("-")

    species_a = species_a.strip()
    species_b = species_b.strip()

    species_tuple = SpeciesTuple(species_a, species_b)

    tokens = value.split()
    if not tokens:
      raise ConfigParserException("[Pair] potential parameter lines should be of the form 'SPECIES_A-SPECIES_B : POTENTIAL_FORM PARAMS...'")

    potential_form = tokens[0]
    params = [float(v) for v in tokens[1:]]

    return PairPotentialTuple(species_tuple, potential_form, params)

  def _parse_potential_form_signature(self, pf):
    pf = pf.strip()
    m = self._signature_re.match(pf)

    if not m:
      raise ConfigParserException("Invalid function signature found in [Potential-Form]: '{0}'".format(pf))

    label, params = m.groups()
    label = label.strip()

    params = [p.strip() for p in params.split(',')]
    return PotentialFormSignatureTuple(label, params)

  @property
  def pair(self):
    """Returns the contents of the config file's [Pair] section.

    :returns: List of tuples of (SpeciesPair, potential_form_label, params)
      Where params = [p1, p2, ..., pn] and p1 etc are the potential parameters"""
    if not self._config_parser.has_section("Pair"):
      raise ConfigParserException("Configuration file does not contain [Pair] section.")
    
    pairparams = []
    for k in self._config_parser['Pair']:
      v = self._config_parser['Pair'][k]
      pair_tuple = self._parse_pair_line(k,v)
      pairparams.append(pair_tuple)
    return pairparams


  @property
  def potential_form(self):
    """Return the contents of the config file's [Potential-Form] section.

    :returns: List of (PotentialFormSignature, formula_string) pairs."""

    if not self._config_parser.has_section("Potential-Form"):
      raise ConfigParserException("Configuration file does not contain [Potential-Form] section.")

    potential_forms = []
    for k in self._config_parser['Potential-Form']:
      v = self._config_parser['Potential-Form'][k]
      signature = self._parse_potential_form_signature(k)
      pf = PotentialFormTuple(signature, v.strip())
      potential_forms.append(pf)
    return potential_forms