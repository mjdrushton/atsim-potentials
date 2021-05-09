import logging

from .._potential import Potential

from ._potential_form_builder import Potential_Form_Builder
from ._potential_form_builder import UnknownModifierException, UnknownPotentialFormException

from ._common import ConfigurationException, Unknown_Modifier_Exception

class Pair_Potentials_From_Tuples_Builder(object):
  """Converts PairPotentialTuple instances into PairPotential objects"""

  def __init__(self, potential_tuples, potential_form_registry, modifier_registry, log_section_name = "Pair"):
    """Create object which builds pair potential instances from tuples.

    Args:
        potential_tuples (list): List of PairPotentialTuple instances
        potential_form_registry (Potential_Form_Registry): Available potential forms
        modifier_registry (Modifier_Registry): Available potential modifiers
        log_section_name (str, optional): Configuration file section name used for logging. Defaults to "Pair".
    """
    self.potential_tuples = potential_tuples
    self.potential_form_registry = potential_form_registry
    self.modifier_registry = modifier_registry
    self.log_section_name = log_section_name
    self._potlist = None

  def _create_potential(self, potrow, mrpfb):
    logger = logging.getLogger(__name__).getChild("Pair_Potentials_From_Tuples_Builder._create_potential")
    logger.debug("Creating potential object for potential row from [{}]: {}".format(self.log_section_name, potrow))
    pot_func = mrpfb.create_potential_function(potrow.potential_form_instance)

    # Make the potential object
    potobj = Potential(potrow.species.species_a, potrow.species.species_b, pot_func)
    return potobj

  def _init_potentials(self):
    pots = []
    pfb = Potential_Form_Builder(self.potential_form_registry, self.modifier_registry)

    for potrow in self.potential_tuples:
      try:
        pot = self._create_potential(potrow, pfb)
        pots.append(pot)
      except UnknownModifierException as ume:
        msg = "Unknown modifier '{modifier_name}' for {species_a}-{species_b} in [{section_name}] section".format(
          modifier_name = ume.args[0],
          species_a = potrow.species.species_a,
          species_b = potrow.species.species_b,
          section_name = self.log_section_name)
        raise Unknown_Modifier_Exception(msg)
      except UnknownPotentialFormException as upe:
        msg = "Unknown potential form '{potform_name}' for {species_a}-{species_b} in [{section_name}] section".format(
          potform_name = upe.args[0],
          species_a = potrow.species.species_a,
          species_b = potrow.species.species_b,
          section_name = self.log_section_name)
        raise ConfigurationException(msg)
      except ConfigurationException as ce:
        msg = "Problem defining {species_a}-{species_b} in [{section_name}] section. {msg}".format(
          msg = ce,
          species_a = potrow.species.species_a,
          species_b = potrow.species.species_b,
          section_name = self.log_section_name)
        raise ConfigurationException(msg)
    return pots

  @property
  def potentials(self):
    if self._potlist is None:
      self._potlist = self._init_potentials()
    return self._potlist

class Pair_Potential_Builder(object):
  """Uses the output of ConfigParser.pair and .potentialforms properties to build
  Potential objects"""

  def __init__(self, cp, potential_form_registry, modifier_registry):
    """:param cp: atsim.potentials.config.ConfigParser instance.
       :param potential_form_register: Potential_Form_Registry
       :param modifier_register: Modifier_Registry"""
    self._tuple_pot_builder = Pair_Potentials_From_Tuples_Builder(cp.pair, potential_form_registry, modifier_registry, "Pair")
    
  @property
  def potentials(self):
    return self._tuple_pot_builder.potentials