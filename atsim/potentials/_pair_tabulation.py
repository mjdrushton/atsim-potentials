from ._lammps_writeTABLE import writePotentials

class LAMMPS_PairTabulation(object):
  """Class for tabulating pair-potential models for LAMMPS"""
  
  def __init__(self, potentials, cutoff, nr):
    """Create pair tabulation for LAMMPS.

    :params potentials: List of atsim.potentials.Potential objects.
    :params cutoff: Maximum separation to be tabulated.
    :params nr: Number of points to be used in tabulation"""
    self._nr = nr
    self._cutoff = cutoff
    self._potentials = potentials
  
  @property
  def type(self):
    return "Pair"

  @property
  def target(self):
    return "LAMMPS"

  @property
  def nr(self):
    return self._nr

  @property
  def cutoff(self):
    return self._cutoff

  @property
  def potentials(self):
    return self._potentials

  @property
  def dr(self):
    return (self.cutoff / float(self.nr-1)) 

  def write(self, fp):
    """Write the tabulation to the file object `fp`.

    :param fp: File object into which data should be written."""
    writePotentials(self.potentials, self.dr, self.cutoff, self.nr-1, fp)
