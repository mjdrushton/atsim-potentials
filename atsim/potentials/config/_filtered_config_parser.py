from ._config_parser import ConfigParser

from wrapt import ObjectProxy


class FilteredConfigParser(ObjectProxy):
  """Class that wraps around ConfigParser instances and
  filters out entries for particular, unwanted species"""


  def __init__(self, config_parser, exclude = [], include = []):
    """Wrap existing ConfigParser so that it excludes entries
    for unwanted species.

    Note: only one of `exclude` and `include` arguments can be specified.

    :param config_parser: ConfigParser instance to be wrapped.
    :param exclude: Collection of species labels, 
      entries for species included in `exclude` will be removed from
      the lists returned by properties such as `pair`, `eam_density` and `eam_embed`.
    :param include: Species labels that should be returned by `pair`, `eam_density` and `eam_embed` functions."""
    ObjectProxy.__init__(self, config_parser)

    if exclude and include:
      raise ValueError("Both exclude and include arguments specified. Only one can be used at one time.")

    if exclude:
      self._species_list = exclude
      self._exclude_flag = True
    else:
      self._species_list = include
      self._exclude_flag = False
    
  def _check_tuple(self, check_tuple):
    for v in check_tuple:
      v_in = v in self._species_list
      if self._exclude_flag and v_in:
        return False
      elif not self._exclude_flag and not v_in:
        return False
    return True

  @property
  def pair(self):
    filtered = [ p for p in self.__wrapped__.pair if self._check_tuple(p.species)]
    return filtered

  @property
  def eam_embed(self):
    filtered = [ p for p in self.__wrapped__.eam_embed if self._check_tuple((p.species,)) ]
    return filtered

  @property
  def eam_density(self):
    filtered = [ p for p in self.__wrapped__.eam_density if self._check_tuple((p.species,)) ]
    return filtered

  @property
  def eam_density_fs(self):
    filtered = [ p for p in self.__wrapped__.eam_density_fs if self._check_tuple(p.species)]
    return filtered
