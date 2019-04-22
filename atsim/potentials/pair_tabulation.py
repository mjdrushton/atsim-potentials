from ._lammps_writeTABLE import writePotentials as lmp_writePotentials
from ._dlpoly_writeTABLE import writePotentials as dlpoly_writePotentials

class _PairTabulation_AbstractBase(object):
  """Base class for PairTabulation objects. 

  Child classes must implement: 
    write() method
    """

  def __init__(self, potentials, cutoff, nr, target):
    """Create pair tabulation for LAMMPS.

    :params potentials: List of atsim.potentials.Potential objects.
    :params cutoff: Maximum separation to be tabulated.
    :params nr: Number of points to be used in tabulation.
    :params target: Label identifying the code for which this object creates tables."""
    self._nr = nr
    self._cutoff = cutoff
    self._potentials = potentials
    self._target = target
  
  @property
  def type(self):
    return "Pair"

  @property
  def target(self):
    return self._target

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
    raise NotImplementedError("Sub-classes must implement write method")


class LAMMPS_PairTabulation(_PairTabulation_AbstractBase):
  """Class for tabulating pair-potential models for LAMMPS"""
  
  def __init__(self, potentials, cutoff, nr):
    """Create pair tabulation for LAMMPS.

    :params potentials: List of atsim.potentials.Potential objects.
    :params cutoff: Maximum separation to be tabulated.
    :params nr: Number of points to be used in tabulation"""
    super(LAMMPS_PairTabulation, self).__init__(potentials, cutoff, nr, "LAMMPS")

  def write(self, fp):
    """Write the tabulation to the file object `fp`.

    :param fp: File object into which data should be written."""
    lmp_writePotentials(self.potentials, self.dr, self.cutoff, self.nr-1, fp)


class DLPoly_PairTabulation(_PairTabulation_AbstractBase):
  """Class for tabulating pair-potential models for DLPOLY"""

  def __init__(self, potentials, cutoff, nr):
    """Create pair tabulation for DLPOLY.

    :params potentials: List of atsim.potentials.Potential objects.
    :params cutoff: Maximum separation to be tabulated.
    :params nr: Number of points to be used in tabulation"""
    super(DLPoly_PairTabulation, self).__init__(potentials, cutoff, nr, "DLPOLY")

  def write(self, fp):
    """Write tabulation to the file object `fp`.

    :param fp: File object into which data should be written."""
    dlpoly_writePotentials(self.potentials, self.cutoff, self.nr, fp)


class GULP_PairTabulation(_PairTabulation_AbstractBase):
  """Class for tabulating pair-potential models for the GULP code.

  .. :seealso::

      * `Gulp Web site <https://nanochemistry.curtin.edu.au/gulp/>`_

  """

  def __init__(self, potentials, cutoff, nr):
    """Create pair tabulation for GULP.

    :params potentials: List of atsim.potentials.Potential objects.
    :params cutoff: Maximum separation to be tabulated.
    :params nr: Number of points to be used in tabulation"""
    super(GULP_PairTabulation, self).__init__(potentials, cutoff, nr, "GULP")

  def write(self, fp):
    """Write tabulation to the file object `fp`.

    :param fp: File object into which data should be written."""
    
    for pot in self.potentials:
      self._write_pot(pot, fp)

  def _write_pot(self, pot, fp):
    header_template = u"{speciesA} {speciesB} {cutoff}\n"
    row_template = u"{energy:.10f} {sepn:.10f}\n"

    fp.write(u"spline cubic\n")
    fp.write(header_template.format(speciesA = pot.speciesA, speciesB = pot.speciesB, cutoff= self.cutoff))

    for n in range(self.nr+1):
      r = float(n)* self.cutoff / (float(self.nr) -1)
      energy = pot.energy(r)
      fp.write(row_template.format(sepn = r, energy = energy))
    fp.write(u"\n")

  