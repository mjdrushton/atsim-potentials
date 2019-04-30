import logging

from .._potential import Potential

from ._potential_form_builder import Potential_Form_Builder
from ._potential_form_builder import UnknownModifierException, UnknownPotentialFormException

from ._common import ConfigurationException, Unknown_Modifier_Exception

class Pair_Potential_Builder(object):
  """Uses the output of ConfigParser.pair and .potentialforms properties to build
  Potential objects"""

  def __init__(self, cp, potential_form_registry, modifier_registry):
    """:param cp: atsim.potentials.config.ConfigParser instance.
       :param potential_form_register: Potential_Form_Registry
       :param modifier_register: Modifier_Registry"""
    self._potlist = self._init_potentials(cp, potential_form_registry, modifier_registry)

  def _init_potentials(self, cp, pfr, mr):
    pots = []
    pfb = Potential_Form_Builder(pfr, mr)

    for potrow in cp.pair:
      try:
        pot = self._create_potential(potrow, pfb)
        pots.append(pot)
      except UnknownModifierException as ume:
        msg = "Unknown modifier '{modifier_name}' for {species_a}-{species_b} in [Pair] section".format(
          modifier_name = ume.args[0],
          species_a = potrow.species.species_a,
          species_b = potrow.species.species_b)
        raise Unknown_Modifier_Exception(msg)
      except UnknownPotentialFormException as upe:
        msg = "Unknown potential form '{potform_name}' for {species_a}-{species_b} in [Pair] section".format(
          potform_name = upe.args[0],
          species_a = potrow.species.species_a,
          species_b = potrow.species.species_b)
        raise ConfigurationException(msg)
      except ConfigurationException as ce:
        msg = "Problem defining {species_a}-{species_b} in [Pair] section. {msg}".format(
          msg = ce,
          species_a = potrow.species.species_a,
          species_b = potrow.species.species_b)
        raise ConfigurationException(msg)
    return pots

  def _create_potential(self, potrow, mrpfb):
    logger = logging.getLogger(__name__).getChild("Pair_Potential_Builder._create_potential")
    logger.debug("Creating potential object for potential row: {}".format(potrow))
    pot_func = mrpfb.create_potential_function(potrow.potential_form_instance)

    # Make the potential object
    potobj = Potential(potrow.species.species_a, potrow.species.species_b, pot_func)
    return potobj

  @property
  def potentials(self):
    return self._potlist