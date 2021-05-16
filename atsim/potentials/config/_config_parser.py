import configparser
import re
import collections

import pyparsing

from ._common import PotentialFormSignatureTuple, \
  PotentialFormTuple, \
  SpeciesTuple, \
  PairPotentialTuple, \
  PotentialFormInstanceTuple, \
  MultiRangeDefinitionTuple, \
  PotentialModifierTuple, \
  EAMFSDensitySpeciesTuple, \
  EAMEmbedTuple, \
  EAMDensityTuple, \
  TableFormTuple

from ._common import ConfigParserException
from ._common import ConfigParserMissingSectionException, ConfigParserDuplicateEntryException

from ._multi_range_parser import multi_range_parser

def _get_or_none(k, d, t):
  v = d.get(k, None)
  if not v is None:
    try:
      v = t(v)
    except ValueError:
      msg = "Could not convert configuration option [{section_name}].{attr_name} into '{type}'. Value is = {value}".format(
        type=t.__name__, 
        section_name = d.name, 
        value = v, 
        attr_name = k)
      raise ConfigParserException(msg)
  return v

ConfigParserOverrideTuple = collections.namedtuple("ConfigParserOverrideTuple", ["section", "key", "value"])

class ConfigOverrideException(ConfigParserException):
  pass

class ConfigOverrideDuplicateException(ConfigOverrideException, ConfigParserDuplicateEntryException):
  pass

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
      nr = int(nr)
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
    'lammps_eam_alloy' : 'setfl',
    'DL_POLY' : 'DLPOLY'
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

  def __repr__(self):
    return "_TabulationSection(r={target}, cutoff={cutoff}, nr={nr}, cutoff_rho={cutoff_rho}, nrho={nrho})".format(
      target = self.target,
      cutoff = self.cutoff,
      nr = self.nr,
      cutoff_rho = self.cutoff_rho,
      nrho = self.nrho)

class _TableFormSection(object):
  """Extracts the [Table-Form:NAME] sections from configuration file.

  These are parsed into TableFormTuple objects which are accessed through the
  `table_forms` property."""

  _section_name_prefix = "Table-Form"
  _section_name_regex = re.compile("^{}:(.*)$".format(_section_name_prefix))

  def __init__(self, cfg_parser):
    self._cfg_parser = cfg_parser
    self._table_forms = self._parse_table_forms()

  def _parse_table_forms(self):
    table_form_list = []
    for section_name in self._cfg_parser.sections():
      if self.is_relevant_section(section_name):
        table_form = self._parse_section(section_name)
        table_form_list.append(table_form)
    return table_form_list

  @classmethod
  def _parse_name(cls, section_name):
    m = cls._section_name_regex.match(section_name)
    name = m.groups()[0]
    name = name.strip()
    return name

  @classmethod
  def is_relevant_section(cls, section_name):
    return cls._section_name_regex.match(section_name) != None

  def _parse_x_y(self, section_name, section):
    x_string = section["x"]
    y_string = section["y"]

    try:
      x = [float(v) for v in x_string.split()]
    except ValueError as e:
      raise ConfigParserException("Error converting value into a float whilst parsing the 'x' entry of '{}': {}".format(section_name, e.args[0]))

    try:
      y = [float(v) for v in y_string.split()]
    except ValueError as e:
      raise ConfigParserException("Error converting value into a float whilst parsing the 'y' entry of '{}': {}".format(section_name, e.args[0]))

    if len(x) != len(y):
      raise ConfigParserException("The number of data items given in the  'x' and 'y' entries of '{}' do not match ({} != {})".format(section_name, len(x), len(y)))

    return (x,y)

  def _parse_xy(self, section_name, section):
    xy_string = section["xy"]

    try:
      xy = [float(v) for v in xy_string.split()]
    except ValueError as e:
      raise ConfigParserException("Error converting value into a float whilst parsing the 'xy' entry of '{}': {}".format(section_name, e.args[0]))

    if len(xy) % 2 != 0:
      raise ConfigParserException("The number of data items in 'xy' is not even for '{}'. This indicates a different number of 'x' and 'y' items".format(section_name))

    even = True
    x = []
    y = []
    for v in xy:
      if even:
        x.append(v)
      else:
        y.append(v)
      even = not even

    return (x,y)

  @classmethod
  def check_for_duplicate_table_forms(cls, cfg_parser):
    seen = {}
    for section_name in cfg_parser.sections():
      if cls.is_relevant_section(section_name):
        label = cls._parse_name(section_name)
        seen.setdefault(label, []).append(section_name)
    
    for _k, v in seen.items():
      if len(v) > 1:
        msg = "Duplicate '{}' sections found:  {}".format(
          cls._section_name_prefix,
          ",".join(["'"+str(s)+"'" for s in v]))
        raise ConfigParserDuplicateEntryException(msg)


  def _parse_data(self, section_name, section):

    # First check that we have reasonable combinations of properties
    # x and y
    # xy

    if "x" in section or "y" in section:
      if not "x" and "y" in section:
        raise ConfigParserException("Did not find both 'x' and 'y' entries whilst parsing the data for section '{}'".format(section_name))

      if "xy" in section:
        raise ConfigParserException("Data in a {} section can either be given using 'xy' or 'x' and 'y' entries. Not both. For section '{}'".format(self._section_name_prefix, section_name))

      data = self._parse_x_y(section_name,section)
    elif "xy" in section:
      if "x" in section or "y" in section:
        raise ConfigParserException("Data in a {} section can either be given using 'xy' or 'x' and 'y' entries. Not both. For section '{}'".format(self._section_name_prefix, section_name))

      data = self._parse_xy(section_name, section)
    else:
      raise ConfigParserException("Could not parse data from '{}', neither 'xy' or 'x' and 'y' entries found.".format(section_name))

    return data

  def _parse_section(self, section_name):
    name = self._parse_name(section_name)
    section = self._cfg_parser[section_name]

    interpolation = section.get(u"interpolation", u"cubic_spline")
    x,y = self._parse_data(section_name, section)

    table_tuple = TableFormTuple(
      name = name, 
      interpolation = interpolation,
      x = x,
      y = y)

    return table_tuple

  @property
  def table_forms(self):
    return self._table_forms
  

class _ConfigParserDict(collections.OrderedDict):
  """Dictionary class used by `_RawConfigParser`,
  this removes whitespace from keys allowing f(x, y) to be equivalent to f(x,y)"""

  def _key_transform(self, k):
    k = k.strip().replace(' ', '')
    k = k.replace('\t', '')
    return k

  def __setitem__(self, key, value):
    key = self._key_transform(key)
    return super(_ConfigParserDict, self).__setitem__(key, value)

  def __getitem__(self, key):
    key = self._key_transform(key)
    return super(_ConfigParserDict, self).__getitem__(key)

  def __delitem__(self, key):
    key = self._key_transform(key)
    return super(_ConfigParserDict, self).__delitem__(key)


class _RawConfigParser(configparser.RawConfigParser):

  def __init__(self):
    super(_RawConfigParser, self).__init__(dict_type = _ConfigParserDict, default_section = "Variables", interpolation = configparser.ExtendedInterpolation())
    self._sections = collections.OrderedDict()

  def optionxform(self, option):
    option = option.strip()
    return option

class ConfigParser(object):
  """Performs initial stage (tokenizing) of generating a potential model
  suitable for tabulation functions."""

  _signature_re = re.compile(r"^([a-zA-Z]\w*?)\((.*)\)")

  # Map of sections relevant to ConfigParser
  # Keys are section keys as the appear to the _config_parser (_RawConfigParser)
  # Values are attribute names on this class to which those sections relate.
  _section_map = {'Tabulation' : 'tabulation',
                'Pair' : 'pair',
                'EAM-Embed' : 'eam_embed',
                'Potential-Form' : 'potential_form',
                'EAM-Density' : None,
                'Table-Form' : 'table_form'}

  def __init__(self, fp, overrides = [], additional = []):
    """Construct ConfigParser.

    :param fp: file object containing configuration data.
    :param overrides: List of `ConfigParserOverrideTuple` instances defining section, key and values of
                      configuration file entries the should be modified before parsing by this class.
    :param additional: List of `ConfigParserOverrideTuple` instances defining section, key and value of entries that should be
                      created in the configuration file before parsing by this class."""
    self._config_parser = self._init_config_parser(fp, overrides, additional)
    self._check_for_duplicates()

    self._tabulation_section = None
    self._table_form = None

    self._default_range_start = MultiRangeDefinitionTuple(u">", 0.0)

  def _init_config_parser(self, fp, overrides, additional):
    cp = _RawConfigParser()
    # cp.readfp(fp)

    try:    
      cp.read_file(fp)
    except (configparser.DuplicateOptionError, configparser.DuplicateSectionError) as e:
      raise ConfigParserDuplicateEntryException(e.message)

    # Process overrides
    for override in overrides:
      if not cp.has_option(override.section, override.key):
        raise ConfigOverrideException(
          "Entry [{section}]: '{key}' not found in configuration file when processing overrides (value = {value})".format(
          section = override.section, key = override.key, value = override.value))

      if override.value is None:
        # Remove item
        cp.remove_option(override.section, override.key)
        if len(cp[override.section]) == 0:
          cp.remove_section(override.section)
      else:
        cp[override.section][override.key] = override.value

    # Add additional values
    for override in additional:
      if cp.has_option(override.section, override.key):
        raise ConfigOverrideDuplicateException(
          "Entry [{section}]: '{key}' already exists in configuration file whilst adding value = {value}".format(
          section = override.section, key = override.key, value = override.value))

      if not cp.has_section(override.section):
        cp.add_section(override.section)
      cp[override.section][override.key] = override.value

    return cp

  def _check_for_duplicates(self):
    self._check_for_duplicate_pairs()
    self._check_for_duplicate_table_forms()

  def _check_for_duplicate_pairs(self):
    """Check the config parser for duplicate pair entries"""

    if self._config_parser.has_section("Pair"):
      seen = set()
      for k in self._config_parser["Pair"]:
        p = self._pair_species_func(k)
        rev_p = tuple(reversed(list(p)))
        if (p in seen) or (rev_p in seen):
          raise ConfigParserDuplicateEntryException("Multiple entries for the pair {A}-{B} found in [Pair] section.".format(A= p[0], B=p[1]))
        seen.add(p)

  def _check_for_duplicate_table_forms(self):
    _TableFormSection.check_for_duplicate_table_forms(self._config_parser)


  def _descend_potential_modifier(self, modifier_node, sibling_iterator, range_defn):

    # ... extract modifier name
    modifier_label = modifier_node['modifier_label']

    # ... extract potential_forms
    params = []
    modifier_parameters = modifier_node['modifier_parameters']
    for p in modifier_parameters:
      ptuple = self._descend_tree(iter(p))
      params.append(ptuple)

    n = self._descend_tree(sibling_iterator)

    ret_tuple = PotentialModifierTuple(modifier = modifier_label, 
      potential_forms = params, 
      start = range_defn, 
      next = n)

    return ret_tuple

  def _descend_potential_description(self, curr_node, sibling_iterator, range_defn):
    # ... extract potential_form
    potential = curr_node['potential_label']
    # ... extract parameters
    parameters = [p for p in curr_node['potential_parameters']]
    # ... are there any more (next)
    n = self._descend_tree(sibling_iterator)
    # ... build tuple
    ret_tuple = PotentialFormInstanceTuple(
      potential_form = potential,
      parameters = parameters,
      start = range_defn,
      next = n)
    return ret_tuple

  def _descend_tree(self, tree_it):
    try:
      first = next(tree_it)
    except StopIteration:
      return None

    name = first.getName()

    # Parse tree can start with a range definition
    # or go straight to a modified instance or potential definition.
    if name == 'range_start':
      range_defn = MultiRangeDefinitionTuple(
        range_type = first['range_type'],
        start = first['start'])
      pot = next(tree_it)
    else:
      range_defn = self._default_range_start
      pot = first

    # pot can be a modified instance or potential definition

    # If pot is 'modified' then create PotentialModifierTuple
    if pot.getName() == 'modifier':
      return self._descend_potential_modifier(pot, tree_it, range_defn)

    # If pot is 'potential_description' then create PotentialFormInstanceTuple
    if pot.getName() == 'potential_description':
      return self._descend_potential_description(pot, tree_it, range_defn)
    
    # Shouldn't get here
    raise Exception("Unknown node type when parsing potential instance description: {}".format(pot.getName()))

  def _parse_multi_range(self, k, value, tuple_type = PairPotentialTuple):
    species = k
    value = value.strip()

    try:
      parse_tree = multi_range_parser.parseString(value, parseAll = True)
    except pyparsing.ParseException as exc:
      msg = "Error when defining potential instance: '{value}' for species = {species}. {msg} (at char: {loc}).".format(
        species = species,
        value = value,
        loc = exc.loc,
        msg = exc.msg)
      raise ConfigParserException(msg)

    tree_it = iter(parse_tree[0])
    potform_instance = self._descend_tree(tree_it)
    return tuple_type(species, potform_instance)

  def _parse_label_type_params_line(self, k, value, parse_key_func, tuple_type):
    if parse_key_func:
      species_tuple = parse_key_func(k)
    else:
      species_tuple = k.strip()
    return self._parse_multi_range(species_tuple, value, tuple_type)

  def _pair_species_func(self, k):
    species_a, species_b = k.split("-")
    species_a = species_a.strip()
    species_b = species_b.strip()
    return  SpeciesTuple(species_a, species_b)


  def _parse_pair_line(self, k, value):
    try:
      return self._parse_label_type_params_line(k, value, self._pair_species_func, PairPotentialTuple)
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

  def _parse_eam_fs_density_line(self, k, value):
    def species_func(k):
      from_species, to_species = k.split("->")
      from_species = from_species.strip()
      to_species = to_species.strip()
      return  EAMFSDensitySpeciesTuple(from_species, to_species)

    try:
      return self._parse_label_type_params_line(k, value, species_func, EAMDensityTuple)
    except ConfigParserException as e:
      section_name = 'EAM-Density'
      raise ConfigParserException("Error in format of Finnis-Sinclair [{section_name}] parameter line: {exc} (should be of the form 'FROM_SPECIES->TO_SPECIES : POTENTIAL_FORM PARAMS...')".format(exc = str(e), section_name = section_name))

  def _parse_potential_form_signature(self, pf):
    pf = pf.strip()
    m = self._signature_re.match(pf)

    if not m:
      raise ConfigParserException("Invalid function signature found in [Potential-Form]: '{0}'".format(pf))

    label, params = m.groups()
    label = label.strip()

    params = [p.strip() for p in params.split(',')]
    return PotentialFormSignatureTuple(label, params, False)

  def _parse_params_section(self, section_name, parse_line_func):
    if not self._config_parser.has_section(section_name):
      raise ConfigParserException("Configuration file does not contain [{section_name}] section.".format(section_name = section_name))
    params = []
    for k in self._config_parser[section_name]:
      v = self._config_parser[section_name][k]
      pair_tuple = parse_line_func(k,v)
      params.append(pair_tuple)
    return params

  def parse_pair_like(self, section_name):
    """Parse a section as if it contains pair potentials.
    
    :param section_name: Name of section that should be parsed in the same way as the [Pair] section.

    :returns: List of tuples of (SpeciesPair, potential_form_label, params)
      Where params = [p1, p2, ..., pn] and p1 etc are the potential parameters"""
    return self._parse_params_section(section_name, self._parse_pair_line)

  @property
  def pair(self):
    """Returns the contents of the config file's [Pair] section.

    :returns: List of tuples of (SpeciesPair, potential_form_label, params)
      Where params = [p1, p2, ..., pn] and p1 etc are the potential parameters"""
    return self.parse_pair_like("Pair")

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
  def table_form(self):
    """Returns parsed content of config file's [Table-Form] section.

    This allows pre-tabulated data to be used within atsim.potentials.

    :returns: List of TableFormTuple instance tuples."""
    if self._table_form is None:
      self._table_form = _TableFormSection(self._config_parser).table_forms
    return self._table_form


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

  @property
  def eam_density_fs(self):
    """Return the parsed contents of the configuration file's [EAM-Density] section.

    This assumes Finnis-Sinclair parsing rules. This means that SPECIES (below)
    is parsed as a `EAMFSDensitySpeciesTuple` with  `from_species` and `to_species` attributes.

    :returns: List of (SPECIES, potential_form_label, params)
      Where params = [p1, p2, ..., pn] and p1 etc are the density function parameters )"""
    return self._parse_params_section("EAM-Density", self._parse_eam_fs_density_line)

  @property
  def parsed_sections(self):
    """Returns a list of relevant sections found inside configuration file. 

    Names are returned as the `ConfigParser` attribute names which could be used to access each parsed section.
    So `[Pair]` becomes `pair` and `[EAM-Density]` is `eam_density`.

    :returns: List of attribute names representing parseable sections of the configuration file"""
    sections = []
    for section_key, output_key in self._section_map.items():
      if output_key and self._config_parser.has_section(section_key):
        sections.append(output_key)

    # The EAM-Density section is overloaded for Finnis-Sinclair and regular EAM potentials
    # if the section's keys contain -> then it's FS
    if self._config_parser.has_section("EAM-Density"):
      isFS = False
      for k in self._config_parser["EAM-Density"]:
        if "->" in k:
          isFS = True
          break
      if isFS:
        sections.append("eam_density_fs")
      else:
        sections.append("eam_density")
    return sections

  @property
  def orphan_sections(self):
    """Returns list of section keys, in current configuration file, that are not relevant to the `ConfigParser` class.property

    :returns: List of section labels."""
    sections = []
    for section_key in self._config_parser.sections():
      if not _TableFormSection.is_relevant_section(section_key ) and not section_key in self._section_map:
        sections.append(section_key)
    return sections

  @property
  def raw_config_parser(self):
    return self._config_parser

  def _convert_species_type(self, property_name, v):
    def default(v):
      return u"{}".format(v)

    known_properties = {
      'atomic_mass' : float, 
      'atomic_number' : int, 
      'covalent_radius' : float,
      'lattice_constant' : float, 
      'charge' : float,
      'lattice_type' : default}

    converted = known_properties.get(property_name, default)(v)
    return converted

  @property
  def species(self):
    """Return reference data for atomic species.
    
    Data is returned as a dictionary relating each species label to a dictionary mapping property name to propety value.
    
    :returns: Dictionary of dictionaries."""
    if not self._config_parser.has_section("Species"):
      return {}

    d = {}
    for k in self._config_parser["Species"]:
      # Split k into SPECIES.PROPERTY_NAME pairs
      tokens = k.split(".", 1)
      if len(tokens) == 1:
        raise ConfigParserException("Error when parsing [Species] section. Keys should be of the form 'SPECIES_LABEL.PROPERTY_NAME'. Invalid key found: '{}'".format(k))
      species, property_name = [t.strip() for t in tokens]
      v = self._config_parser["Species"][k]
      v = self._convert_species_type(property_name, v)
      d.setdefault(species, {})[property_name] = v
    return d