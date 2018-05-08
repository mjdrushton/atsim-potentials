from backports import configparser
import re


from ._common import PotentialFormSignatureTuple, PotentialFormTuple, SpeciesTuple, PairPotentialTuple
from ._common import ConfigParserException
from ._common import ConfigParserMissingSectionException


class _TabulationSection(object):
  """Represents the [Tabulation] section of a config file"""

  def __init__(self, cp):
    self._target = None
    self._type = None
    self._cutoff = None
    self._nr = None
    self._init_from_config(cp)

  def _init_from_config(self, cp):
    # Tabulation section not defined
    if not cp.has_section('Tabulation'):
      return

    self._init_cutoff(cp)

  def _get_or_none(self, k, d, t):
    v = d.get(k, None)
    if not v is None:
      v = t(v)
    return v

  def _init_cutoff(self, cp):
    nr = self._get_or_none('nr', cp['Tabulation'], int)
    dr = self._get_or_none('dr', cp['Tabulation'], float)
    cutoff = self._get_or_none('cutoff', cp['Tabulation'], float)

    if nr and dr and cutoff:
      raise ConfigParserException("'cutoff', 'nr' and 'dr' cannot all be spcified in [Tabulation] section of potential definition.")
    elif nr and dr:
      # Set cutoff
      cutoff = (nr-1)*dr      
    elif cutoff and dr:
      # Set nr
      nr = (cutoff/dr) + 1
    elif not dr is None:
      raise ConfigParserException("'dr' cannot be specified without either 'nr' or 'cutoff' in [Tabulation] section of potential definition.")

    if not nr is None and nr <= 0:
      raise ConfigParserException("'nr' in [Tabulation] section of potential definition cannot be 0 (zero) or negative.")
    if not dr is None and dr <= 0:
      raise ConfigParserException("'dr' in [Tabulation] section of potential definition cannot be 0 (zero) or negative.")
    if not cutoff is None and cutoff <= 0:
      raise ConfigParserException("'cutoff' in [Tabulation] section of potential definition cannot be 0 (zero) or negative.")

    self._nr = nr 
    self._cutoff = cutoff

  @property
  def target(self):
    return self._target

  @property
  def type(self):
    return self._type

  @property
  def cutoff(self):
    return self._cutoff

  @property
  def nr(self):
    return self._nr

class _RawConfigParser(configparser.RawConfigParser):

  def optionxform(self, option):
    return option

class ConfigParser(object):
  """Performs initial stage (tokenizing) of generating a potential model
  suitable for tabulation functions."""

  _signature_re = re.compile(r"^([a-zA-Z]\w*?)\((.*)\)")

  def __init__(self, fp):
    """Construct ConfigParser.

    :param fp: file object containing configuration data."""
    self._config_parser = self._init_config_parser(fp)
    self._tabulation_section = None

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
      raise ConfigParserMissingSectionException("Configuration file does not contain [Potential-Form] section.")

    potential_forms = []
    for k in self._config_parser['Potential-Form']:
      v = self._config_parser['Potential-Form'][k]
      signature = self._parse_potential_form_signature(k)
      pf = PotentialFormTuple(signature, v.strip())
      potential_forms.append(pf)
    return potential_forms

  @property
  def tabulation(self):
    """Return the parsed contents of the config file's [Tabulation] section.

    This defines what type of model (pair, EAM) the config file contains and also how
    the model should be tabulated.

    :returns: _Tabulation_Section object"""
    if self._tabulation_section is None:
      self._tabulation_section = _TabulationSection(self._config_parser)
    return self._tabulation_section