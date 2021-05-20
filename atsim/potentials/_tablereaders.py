import re
import bisect

class _XProxy(object):
  """Class that proxies list of (x,y) tuples returning only the x component.
  Used by TableReaderBase to allow use of the bisect module"""

  def __init__(self, wrapped):
    self._wrapped = wrapped

  #def __getslice__(self, slice):
    #for x,y in self._wrapped[slice]:
    #yield x

  def __getitem__(self, idx):
    return self._wrapped[idx][0]

  def __len__(self):
    return len(self._wrapped)

class TableReaderBase(list):
  """Abstract base class for linearly interpolated base classes, classes inheriting
  from TableReaderBase must implement _populate method that accepts a file object
  and populates self (TableReaderBase inherits from list), with ordered (x,y) tuples"""

  def __init__(self, fileobj, inputConvert = None, outputConvert = None):
    self._populate(fileobj)
    self.xproxy = _XProxy(self)

    if inputConvert == None and outputConvert == None:
      return

    if inputConvert == None:
      inputConvert = lambda x: x

    if outputConvert == None:
      outputConvert = lambda x: x

    for (i, (x,y)) in enumerate(self):
      self[i] = (inputConvert(x), outputConvert(y))


  def _populate(self, fileobj):
    raise Exception("Not implemented")

  def getValue(self, x):
    """Returns y value for given x value, using linear interpolation
    if x lies between data-points. If x is out of the range described
    by this object then returns 0.0

    @param x x-value for which y value should be calculated
    @return y value for given x-position"""
    lowidx = self._findIndex(x)
    if lowidx == None:
      return 0.0

    lx, ly = self[lowidx]
    if lx == x:
      return ly

    highidx = lowidx +1
    if highidx == len(self):
      return 0.0

    hx,hy = self[highidx]

    m = (hy-ly)/(hx - lx)
    c = ly - (m*lx)
    return (m*x) + c

  def _findIndex(self, x):
    """Returns the index of the last x value in this object that is less than x.

    @param x Value for which index should be found
    @return Index into this collection for list value that has x component less than search x"""
    if x< self[0][0] or x> self[-1][0]:
      return None

    idx = bisect.bisect_left(self.xproxy, x)
    if self[idx][0] == x:
      return idx
    else:
      return idx-1


class DatReader(TableReaderBase):
  """Class that reads space separated data from a text file providing linearly interpolated
  value lookups through getValue() method"""

  def _populate(self, fileobj):
    splitre = re.compile(r'\s+')
    results = []
    for line in fileobj:
      line = line[:-1]
      line = line.strip()
      if len(line) == 0 or line[0] == '#':
        continue
      (x,y) = splitre.split(line)[:2]

      results.append( (float(x), float(y) ))
    results.sort()
    self.extend(results)
