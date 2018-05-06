import collections
from backports import configparser
import re
import itertools

import cexprtk

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


class Potential_Form_Registry_Exception(Exception):
  pass

class Potential_Form_Registry(object):
  """Factory class that takes [Potential-Form] definitions
  from ConfigParser and turns them into Potential_Form objects"""

  def __init__(self, cfg):
    definitions = cfg.potential_form
    self._potential_forms = self._build_potential_forms(definitions)
    self._register_with_each_other()
    self._definitions = definitions


  def _build_potential_forms(self, definitions):
    potential_forms = {}
    for d in definitions:
      if d.signature.label in potential_forms:
        raise Potential_Form_Registry_Exception("Two potential forms have the same label: '{0}'".format(d.signature.label))
      pf = Potential_Form(d)
      potential_forms[d.signature.label] = pf
    return potential_forms

  def _register_with_each_other(self):
    # So that each function can rely on other custom functions, add each function to every other
    # function's symbol table
    pairs = list(itertools.permutations(self._potential_forms.values(), 2))
    # print([(a.signature.label, b.signature.label) for (a,b) in pairs])
    for a,b in pairs:
      a._function.register_function(b._function)

  @property
  def registered(self):
    """Returns the labels for the potentials registered here."""
    return [d.signature.label for d in self._definitions]

  def __getitem__(self, k):
    return self._potential_forms[k]

class Potential_Form_Exception(Exception):
  pass

class _Potential_Function(object):
  """Callable that can be added to cexprtk symbol_table. 

  Its wrapped cexprtk expression is only instantiated on first use. This is to allow
  the expression's symbol table to be populated before the expression is parsed. 
  Otherwise custom functions may not have been defined at the time of parsing."""


  def __init__(self, potential_form_tuple):
    """:param potential_form_tuple: PotentialFormTuple describing this function"""
    self._potential_form_tuple = potential_form_tuple
    self._local_symbol_table = self._init_symbol_table()
    self._expression = None

  def _init_symbol_table(self):
    local_symbol_table = cexprtk.Symbol_Table({})
    parameter_names = self._potential_form_tuple.signature.parameter_names
    for pn in parameter_names:
      local_symbol_table.variables[pn] = 1.0
    return local_symbol_table

  def register_function(self, func):
    """Register `func` with this object's symbol_table"""
    label = func._potential_form_tuple.signature.label
    self._local_symbol_table.functions[label] = func

  def __call__(self, *args):
    parameter_names = self._potential_form_tuple.signature.parameter_names
    assert len(args) == len(parameter_names)
    for (pn, v) in zip(parameter_names, args):
      self._local_symbol_table.variables[pn] = v

    if not self._expression:
      self._expression = cexprtk.Expression(self._potential_form_tuple.expression, self._local_symbol_table)
    retval = self._expression()
    return retval

class Potential_Form(object):
  
  def __init__(self, potential_form):
    """Create Potential_Form object.

    :param potential_form: PotentialFormTuple describing the potential."""
    self.potential_definition = potential_form
    self._function = _Potential_Function(potential_form)

  @property
  def signature(self):
    return self.potential_definition.signature

  @property
  def expression(self):
    return self.potential_definition.expression

  def __call__(self, *args):
    # Curry all but the first argument
    if not len(args) == len(self.signature.parameter_names)-1:
      raise Potential_Form_Exception(
        "Potential form requires {0} arguments but {1} were provided".format(len(self.signature.parameter_names)-1, len(args)))

    def f(r):
      funcargs = [r]
      funcargs.extend(args)
      return self._function(*funcargs)
    return f