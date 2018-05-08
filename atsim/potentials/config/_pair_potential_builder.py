import logging

from .._potential import Potential

class Pair_Potential_Builder(object):
  """Uses the output of ConfigParser.pair and .potentialforms properties and builds
  Potential objects"""

  def __init__(self, cp, potential_form_registry):
    """:param cp: atsim.potentials.config.ConfigParser instance.
       :param potential_form_register: Potential_Form_Registry"""
    self._potlist = self._init_potentials(cp, potential_form_registry)

  def _init_potentials(self, cp, pfr):
    pots = []
    for potrow in cp.pair:
      pot = self._create_potential(potrow, pfr)
      pots.append(pot)
    return pots

  def _create_potential(self, potrow, pfr):
    logger = logging.getLogger(__name__).getChild("Pair_Potential_Builder._create_potential")
    logger.debug("Creating potential object for potential row: {}".format(potrow))
    potential_form = pfr[potrow.potential_form]
    
    # Parametrise the potential form to create a potential function
    pot_func = potential_form(*potrow.parameters)

    # Make the potential object
    potobj = Potential(potrow.species.species_a, potrow.species.species_b, pot_func)
    return potobj

  @property
  def potentials(self):
    return self._potlist