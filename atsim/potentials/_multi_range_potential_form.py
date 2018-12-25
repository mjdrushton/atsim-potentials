# import collections
import functools
import operator

from ._util import gradient

def _range_defn_cmp(a,b):
  if a.start == b.start:
    if a.range_type == ">=" and b.range_type == ">":
      return -1
    if b.range_type == ">=" and a.range_type == ">":
      return 1
    return 0

  return (a.start > b.start) - (a.start < b.start)

_range_defn_key = functools.cmp_to_key(_range_defn_cmp)

class Multi_Range_Defn(object):

  def __init__(self, range_type, start, potential_form, **kwargs):
    self._range_type = range_type
    self._start = start
    self._potential_form = potential_form

    self._deriv_callable = gradient(self._potential_form)
    self._deriv2_callable = gradient(self._deriv_callable)

  @property
  def range_type(self):
    return self._range_type

  @property
  def start(self):
    return self._start

  @property
  def potential_form(self):
    return self._potential_form

  @property
  def has_deriv(self):
    """Returns True if the potential callable provides an analytical derivative through a `.deriv()` method."""
    return hasattr(self.potential_form, "deriv")

  @property
  def has_deriv2(self):
    """Returns True if the potential callable provides an analytical derivative through a `.deriv2()` method."""
    return hasattr(self.potential_form, "deriv2")

  def deriv(self, r):
    return self._deriv_callable(r)

  def deriv2(self, r):
    return self._deriv2_callable(r)

def create_Multi_Range_Potential_Form(*range_tuples, **kwargs):
  """Creates Multi_Range_Potential_Form or sub-class instance, from list of Multi_Range_Defn
  instances in `range_tuples`.

  If any Multi_Range_Defn object's `.has_deriv2` are True then an instance
  of Multi_Range_Potential_Form_Deriv2 is returned.

  If any Multi_Range_Defn object's `.has_deriv` property is True but all `.has_deriv2` are False then an instance
  of Multi_Range_Potential_Form_Deriv is returned.

  If non of the Multi_Range_Defn objects provide analytical deriv or deriv2 methods, return Multi_Range_Potential_Form.

  :param range_tuples: List of Multi_Range_Defn instances.
  :param kwargs: Keyword arguments passed to Multi_Range_Potential_Form constructor.

  :return: See above"""

  any_deriv2 = False
  any_deriv = False
  for rt in range_tuples:
    any_deriv2 = any_deriv2 or rt.has_deriv2
    any_deriv = any_deriv or rt.has_deriv

  if any_deriv2:
    cls =  Multi_Range_Potential_Form_Deriv2
  elif any_deriv:
    cls =  Multi_Range_Potential_Form_Deriv
  else:
    cls =  Multi_Range_Potential_Form
  
  obj = cls(*range_tuples, **kwargs)

  return obj


class Multi_Range_Potential_Form(object):
  """Class allowing the creation of composite potential forms. 
  This allows different potential functions to be defined for different 
  distance ranges when the potential is called"""

  def __init__(self, *range_defns, **kwargs):
    """Define potential form from a list of `Multi_Range_Defn_Tuple` objects.

    Each named tuple provides a callable (each accepting a single argument) that should be used for a particular range.
    The `start` property of the tuple defines the point at which this callable may be used.
    The `range_type` property describes how the boundary between two ranges should be handled,
    this can be either `>` or `>=`. If `>` then only use that range's callable when *above* the `start` value.
    If `>=` then the range is inclusive of the `start` value.

    Example:

      To define a two range potential form where potential where the callable `func1` is used
      for values less than 2.0 and `func2` beyond this the following definition may be used:

      .. code-block::

          multi_range_potential = Multi_Range_Potential_Form(
            Multi_Range_Defn(">", float("-inf'), func1),
            Multi_Range_Defn(">", 2.0, func2))

    
    :param range_tuples: List of Multi_Range_Defn defining potential ranges.
    :param default_value: This value is returned when this object is called with an argument below the 
      the lowest `start` value in `range_tuples`."""

    if kwargs:
      keys = set(kwargs)
      keys = keys - set(['default_value'])
      if keys:
        raise ValueError("Unknown keyword arguments: {}".format(",".join(sorted(keys))))

    self.default_value = kwargs.get('default_value', 0.0)
    self._range_defns = None
    self.range_defns = range_defns

  def _range_search(self, r):
    rt = self.range_defns
    if not rt or r < rt[0].start or (r == rt[0].start and rt[0].range_type == '>'):
      return None

    last = None
    for t in rt:
      if r == t.start and t.range_type == '>=':
        return t
      elif last and r <= t.start and r > last.start:
        return last
      else:
        last = t
    
    if last and r > last.start:
      return last

  @property
  def range_defns(self):
    return self._range_defns

  @range_defns.setter
  def range_defns(self, range_defns):
    tuples = list(range_defns)
    tuples.sort(key = _range_defn_key )
    self._range_defns = tuples

  def __call__(self, r):
    rt = self._range_search(r)

    if rt is None:
      return self.default_value
    
    return rt.potential_form(r)


class Multi_Range_Potential_Form_Deriv(Multi_Range_Potential_Form):
  """Sub-class of Multi_Range_Potential_Form which additionally provides .deriv() method"""


  def deriv(self, r):
    rt = self._range_search(r)
    if rt is None:
      return 0.0
    return rt.deriv(r)

class Multi_Range_Potential_Form_Deriv2(Multi_Range_Potential_Form_Deriv):
  """Sub-class of Multi_Range_Potential_Form which additionally provides .deriv() and .deriv2() methods"""

  def deriv2(self, r):
    rt = self._range_search(r)
    if rt is None:
      return 0.0
    return rt.deriv2(r)