from ._config_parser import ConfigParser

from ._common import ConfigurationException

from ._tabulation_factories import TABULATION_FACTORIES

import logging

class Configuration(object):
  """Factory class that allows Tabulation objects to be built from .ini files"""

  def __init__(self):
    self._tabulation_factories = dict(TABULATION_FACTORIES)

  def read(self, fp):
    """Read potential data from the file object `fp` and return a `PairTabulation` or `EAMTabulation` object.

    :params fp: File like object containing potential information.
    :returns: Tabulation object"""
    cp = ConfigParser(fp)
    return self.read_from_parser(cp)

  def read_from_parser(self, cp):
    """Read potential data from the `ConfigParser` object `cp` and return a `PairTabulation` or `EAMTabulation` instance.

    :param cp: atsim.potentials.config.ConfigParser instance.
    :returns: Tabulation object"""
    logger = logging.getLogger(__name__).getChild("Configuration.read_from_parser")

    tabulation_target = cp.tabulation.target
    if tabulation_target is None:
      logger.warning("No tabulation target specified - defaulting to 'LAMMPS'")
      tabulation_target = "LAMMPS"

    if not tabulation_target in self._tabulation_factories:
      errormsg = "[Tabulation].tabulation-target - unknown tabulation target specified : '{}'".format(tabulation_target)
      logger.error(errormsg)
      raise ConfigurationException(errormsg)
    else:
      logger.info("Tabulation target specified as '{}'".format(tabulation_target))

    factory = self._tabulation_factories[tabulation_target]
    tabulation = factory.create_tabulation(cp)
    return tabulation
    

    