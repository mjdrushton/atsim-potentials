from .._modifiers import is_modifier

import logging
import inspect


class Modifier_Registry(object):
  """Registry of factories for potential modifiers"""

  def __init__(self):
    self._logger = logging.getLogger(__name__).getChild("Modifier_Registry")
    self._modifiers = self._register_standard()

  def _register_standard(self):
    logger = self._logger.getChild("_register_standard")
    from .. import _modifiers
    modifiers = {}
    for name, pyobj in inspect.getmembers(_modifiers, is_modifier):
      logger.debug("Registering modifier: {}".format(name))
      modifiers[name] = pyobj
    return modifiers
  
  def __getitem__(self, k):
    return self._modifiers[k]