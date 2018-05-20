from ._config_parser import ConfigParser
from ._pair_potential_builder import Pair_Potential_Builder
from ._potential_form_registry import Potential_Form_Registry

from ._common import PAIR_TABULATION

from .._pair_tabulation import LAMMPS_PairTabulation

import logging

class Configuration(object):
  """Factory class that allows Tabulation objects to be built from .ini files"""
  
  def read(self, fp):
    """Read potential data from the file object `fp` and return a `PairTabulation` or `EAMTabulation` object.

    :params fp: File like object containing potential information.
    :returns: Tabulation object"""
    logger = logging.getLogger(__name__).getChild("Configuration.read")
    cp = ConfigParser(fp)

    tabulation_type = cp.tabulation.type
    if tabulation_type is None:
      logger.debug("No tabulation type specified - defaulting to 'Pair'")
      tabulation_type = PAIR_TABULATION

    if tabulation_type == PAIR_TABULATION:
      return self._build_pair_tabulation(cp)
    else:
      raise ValueError("Unknown tabulation type '{}'".format(tabulation_type))

  def _create_pair_objects(self, cp):
    potential_form_registry = Potential_Form_Registry(cp, True)
    potbuilder = Pair_Potential_Builder(cp, potential_form_registry)
    pots = potbuilder.potentials
    return pots

  def _build_pair_tabulation(self, cp):
    logger = logging.getLogger(__name__).getChild("Configuration._build_pair_tabulation")

    tabulation_target = cp.tabulation.target
    if tabulation_target is None:
      logger.warning("No tabulation target specified - defaulting to 'LAMMPS'")
      tabulation_target = "LAMMPS"

    # Get cutoff and gridpoints
    if cp.tabulation.cutoff is None:
      cutoff = 10.0
    else:
      cutoff = cp.tabulation.cutoff
    
    if cp.tabulation.nr is None:
      nr = 1001
    else:
      nr = cp.tabulation.nr

    # Get pair potentials
    potobjs = self._create_pair_objects(cp)

    logger.info("Creating pair-potential tabulation with following properties:")
    logger.info("  * tabulation-target : {}".format(tabulation_target))
    logger.info("  * cutoff: {}".format(cutoff))
    logger.info("  * nr: {}".format(nr))
    logger.info("  * Tabulation will contain:")
    for pot in potobjs:
      logger.info("      + {}-{}".format(pot.speciesA, pot.speciesB))

    if tabulation_target == "LAMMPS":
      # Make LAMMPS tabulation here
      tabulation = LAMMPS_PairTabulation(potobjs, cutoff, nr)
      return tabulation
    else:
      raise ValueError("Unknown pair tabulation target '{}'".format(tabulation_target))

      