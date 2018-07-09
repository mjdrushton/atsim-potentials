import logging

from .._potential import Potential

from ._potential_form_builder import Potential_Form_Builder

class Pair_Potential_Builder(object):
  """Uses the output of ConfigParser.pair and .potentialforms properties and builds
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
      pot = self._create_potential(potrow, pfb)
      pots.append(pot)
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