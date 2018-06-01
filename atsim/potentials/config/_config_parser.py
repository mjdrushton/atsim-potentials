from backports import configparser
import re
import collections

from ._common import PotentialFormSignatureTuple, PotentialFormTuple, SpeciesTuple, PairPotentialTuple, PotentialFormInstanceTuple, MultiRangeDefinitionTuple
from ._common import EAMEmbedTuple, EAMDensityTuple
from ._common import ConfigParserException
from ._common import ConfigParserMissingSectionException

def _get_or_none(k, d, t):
  v = d.get(k, None)
  if not v is None:
    v = t(v)
  return v

class _TabulationCutoff(object):

  def __init__(self, cutoff_name, nr_attr='nr', dr_attr = 'dr', cutoff_attr = 'cutoff'):
    self._cutoff_name = cutoff_name
    self._nr_attr = nr_attr
    self._dr_attr = dr_attr
    self._cutoff_attr = cutoff_attr
    self._template_dict = dict( cutoff = self._cutoff_attr, nr = self._nr_attr, dr = self._dr_attr) 

  def create_cutoff(self, cp_tabulation_section):
    tclass = collections.namedtuple(self._cutoff_name, [self._nr_attr, self._cutoff_attr])
    nr = None
    cutoff = None
    if cp_tabulation_section:
      nr,cutoff = self._init_cutoff(cp_tabulation_section)
    return tclass(nr, cutoff)

  def _init_cutoff(self, cp_tabulation_section):
    nr = _get_or_none(self._nr_attr, cp_tabulation_section, int)
    dr = _get_or_none(self._dr_attr, cp_tabulation_section, float)
    cutoff = _get_or_none(self._cutoff_attr, cp_tabulation_section, float)

    if nr and dr and cutoff:
      raise ConfigParserException("'{cutoff}', '{nr}' and '{dr}' cannot all be spcified in [Tabulation] section of potential definition.".format(**self._template_dict))
    elif nr and dr:
      # Set cutoff
      cutoff = (nr-1)*dr      
    elif cutoff and dr:
      # Set nr
      nr = (cutoff/dr) + 1
    elif not dr is None:
      raise ConfigParserException("'{dr}' cannot be specified without either '{nr}' or '{cutoff}' in [Tabulation] section of potential definition.".format(**self._template_dict))

    if not nr is None and nr <= 0:
      raise ConfigParserException("'{nr}' in [Tabulation] section of potential definition cannot be 0 (zero) or negative.".format(**self._template_dict))
    if not dr is None and dr <= 0:
      raise ConfigParserException("'{dr}' in [Tabulation] section of potential definition cannot be 0 (zero) or negative.".format(**self._template_dict))
    if not cutoff is None and cutoff <= 0:
      raise ConfigParserException("'{cutoff}' in [Tabulation] section of potential definition cannot be 0 (zero) or negative.".format(**self._template_dict))
    return nr, cutoff

class _TabulationSection(object):
  """Represents the [Tabulation] section of a config file"""

  _target_synonyms = {
    'lammps_eam_alloy' : 'setfl'
  }

  def __init__(self, cp):
    self._target = None
    self._cutoff = None
    self._nr = None
    self._init_from_config(cp)

  def _init_from_config(self, cp):
    self._init_cutoff(cp)
    # Tabulation section not defined
    if not cp.has_section('Tabulation'):
      return
    self._init_target(cp)

  def _init_cutoff(self, cp):
    if cp.has_section('Tabulation'):
      tabulation_section = cp['Tabulation']
    else:
      tabulation_section = None
    r_cutoff = _TabulationCutoff('R_Cutoff').create_cutoff(tabulation_section)
    density_cutoff = _TabulationCutoff('Density_Cutoff', 'nrho', 'drho', 'cutoff_rho').create_cutoff(tabulation_section)

    self._r_cutoff = r_cutoff
    self._density_cutoff = density_cutoff

  def _init_target(self, cp):
    target = _get_or_none('target', cp['Tabulation'], str)

    # For convenience some targets can be specified with more than one value.
    # these are stored in the _target_synonums dict. Check this now.
    target = self._target_synonyms.get(target, target)
    self._target = target

  @property
  def target(self):
    return self._target

  @property
  def cutoff(self):
    return self._r_cutoff.cutoff

  @property
  def nr(self):
    return self._r_cutoff.nr

  @property
  def cutoff_rho(self):
    return self._density_cutoff.cutoff_rho

  @property
  def nrho(self):
    return self._density_cutoff.nrho

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

  def _parse_multi_range(self, k, value, tuple_type = PairPotentialTuple):
    species = k
    regex = "({cap}>=?)\s*({cap}[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)"
    splitter = "({})".format(regex.format(cap = "?:"))
    matcher = regex.format(cap = "")
    
    #TODO: Store the compiled regexes on the class.
    split_regex = re.compile(splitter)
    match_regex = re.compile(matcher)

    # import pdb; pdb.set_trace()
    value = value.strip()
    isrange = False

    parts = []
    range_tuple = MultiRangeDefinitionTuple(">", 0.0)
    split_parts = split_regex.split(value)
    for i,sub in enumerate(split_parts):
      if i == 0 and not sub:
        isrange = True
        continue
      
      if isrange:
        m = match_regex.match(sub)
        if not m:
          raise ConfigParserException("Could not parse multi-range potential form. '{}' could not be parsed as range for potential-form definition '{} : {}'".format(sub, species, value))
        groups = m.groups()
        range_type = groups[0]
        range_value = groups[1]

        try:
          range_value = float(range_value)
        except ValueError:
          raise ConfigParserException("Could not parse multi-range potential form. '{}' could not be parsed as range for potential-form definition '{} : {}'".format(sub, species, value))
        range_tuple = MultiRangeDefinitionTuple(range_type, range_value)
        isrange = False
      else:
        isrange = True
        tokens = sub.split()

        
        potential_form_label = tokens[0]
        try:
          params = [float(t) for t in tokens[1:]]
        except ValueError as e:
          raise ConfigParserException("Could not parse multi-range potential form. Parameter could not be converted to float for part #{rangenum} '{potential_form_label}': {args}".format(
            rangenum = i+1, 
            potential_form_label = potential_form_label, 
            args = e.args))

        potform_tuple = PotentialFormInstanceTuple(potential_form_label, params, range_tuple, None)
        parts.append(potform_tuple)

    # Now assemble into potential form tuples
    if len(parts) == 1:
      potform_instance = parts[0]
    else:
      last = parts.pop()
      while parts:
        curr = parts.pop()
        last = PotentialFormInstanceTuple(curr.potential_form, curr.parameters, curr.start, last)
      potform_instance = last
    
    return tuple_type(species, potform_instance)

  def _parse_label_type_params_line(self, k, value, parse_key_func, tuple_type):
    if parse_key_func:
      species_tuple = parse_key_func(k)
    else:
      species_tuple = k.strip()
    return self._parse_multi_range(species_tuple, value, tuple_type)

  def _parse_pair_line(self, k, value):
    def species_func(k):
      species_a, species_b = k.split("-")
      species_a = species_a.strip()
      species_b = species_b.strip()
      return  SpeciesTuple(species_a, species_b)

    try:
      return self._parse_label_type_params_line(k, value, species_func, PairPotentialTuple)
    except ConfigParserException:
      raise ConfigParserException("[Pair] potential parameter lines should be of the form 'SPECIES_A-SPECIES_B : POTENTIAL_FORM PARAMS...'")

  def _parse_eam_line(self, k, value, section_name, tuple_type):
    try:
      return self._parse_label_type_params_line(k, value, None, tuple_type)
    except ConfigParserException:
      raise ConfigParserException("[{section_name}] parameter lines should be of the form 'SPECIES : POTENTIAL_FORM PARAMS...'".format(section_name = section_name))

  def _parse_embed_line(self, k, value):
    return self._parse_eam_line(k, value, "EAM-Embed", EAMEmbedTuple)

  def _parse_density_line(self, k, value):
    return self._parse_eam_line(k, value, "EAM-Density", EAMDensityTuple)

  def _parse_potential_form_signature(self, pf):
    pf = pf.strip()
    m = self._signature_re.match(pf)

    if not m:
      raise ConfigParserException("Invalid function signature found in [Potential-Form]: '{0}'".format(pf))

    label, params = m.groups()
    label = label.strip()

    params = [p.strip() for p in params.split(',')]
    return PotentialFormSignatureTuple(label, params)

  def _parse_params_section(self, section_name, parse_line_func):
    if not self._config_parser.has_section(section_name):
      raise ConfigParserException("Configuration file does not contain [{section_name}] section.".format(section_name = section_name))
    params = []
    for k in self._config_parser[section_name]:
      v = self._config_parser[section_name][k]
      pair_tuple = parse_line_func(k,v)
      params.append(pair_tuple)
    return params

  @property
  def pair(self):
    """Returns the contents of the config file's [Pair] section.

    :returns: List of tuples of (SpeciesPair, potential_form_label, params)
      Where params = [p1, p2, ..., pn] and p1 etc are the potential parameters"""
    return self._parse_params_section("Pair", self._parse_pair_line)

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

  @property
  def eam_embed(self):
    """Return the parsed contents of the configuration file's [EAM-Embed] section.

    :returns: List of (SPECIES, potential_form_label, params)
      Where params = [p1, p2, ..., pn] and p1 etc are the embedding function parameters )"""
    return self._parse_params_section("EAM-Embed", self._parse_embed_line)

  @property
  def eam_density(self):
    """Return the parsed contents of the configuration file's [EAM-Density] section.

    :returns: List of (SPECIES, potential_form_label, params)
      Where params = [p1, p2, ..., pn] and p1 etc are the density function parameters )"""
    return self._parse_params_section("EAM-Density", self._parse_density_line)
