
import collections
import functools

def _tuple_cmp(a,b):
  if a.start == b.start:
    if a.range_type == ">=" and b.range_type == ">":
      return -1
    if b.range_type == ">=" and a.range_type == ">":
      return 1
    return 0

  return (a.start > b.start) - (a.start < b.start)

_tuple_key = functools.cmp_to_key(_tuple_cmp)


class Multi_Range_Potential_Form(object):
  """Class allowing the creation of composite potential forms. 
  This allows different potential functions to be defined for different 
  distance ranges when the potential is called"""

  Range_Defn_Tuple = collections.namedtuple("Range_Defn_Tuple", ["range_type", "start", "potential_form"])  

  def __init__(self, *range_tuples, **kwargs):
    """Define potential form from a list of `Multi_Range_Potential_Form.Range_Defn_Tuple` objects.

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
            Multi_Range_Potential_Form.Range_Defn_Tuple(">", float("-inf'), func1),
            Multi_Range_Potential_Form.Range_Defn_Tuple(">", 2.0, func2))

    
    :param range_tuples: List of Range_Defn_Tuple defining potential ranges.
    :param default_value: This value is returned when this object is called with an argument below the 
      the lowest `start` value in `range_tuples`."""

    if kwargs:
      keys = set(kwargs)
      keys = keys - set(['default_value'])
      if keys:
        raise ValueError("Unknown keyword arguments: {}".format(",".join(sorted(keys))))

    self.default_value = kwargs.get('default_value', 0.0)
    self._range_tuples = None
    self.range_tuples = range_tuples

  def _range_search(self, r):
    rt = self.range_tuples
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
  def range_tuples(self):
    return self._range_tuples

  @range_tuples.setter
  def range_tuples(self, range_tuples):
    tuples = list(range_tuples)
    tuples.sort(key = _tuple_key )
    self._range_tuples = tuples

  def __call__(self, r):
    rt = self._range_search(r)

    if rt is None:
      return self.default_value

    return rt.potential_form(r)
