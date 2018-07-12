
from ._pair_tabulation import _PairTabulation_AbstractBase

from ._lammpsWriteEAM import writeSetFL, writeSetFLFinnisSinclair
from ._dlpoly_writeTABEAM import writeTABEAM, writeTABEAMFinnisSinclair

class _EAMTabulationAbstractbase(_PairTabulation_AbstractBase):
  """Base class for EAMTabulation objects.

  Child classes must implement:
    write() method.

  """

  def __init__(self, potentials, eam_potentials, cutoff, nr, cutoff_rho, nrho, target):
    """
    :params potentials: List of atsim.potentials.Potential objects.
    :params eam_potentials: List of `atsim.potentials.EAMPotential` instances.
    :params cutoff: Maximum separation to be tabulated.
    :params nr: Number of points to be used in tabulation
    :params cutoff_rho: Density cutoff.
    :params nrho: Number of points to be used when discretising density range during EAM tabulation
    :params target: Name of tabulation target.
    """
    super(_EAMTabulationAbstractbase, self).__init__(potentials, cutoff, nr, target)
    self._nrho = nrho
    self._cutoff_rho = cutoff_rho
    self._eam_potentials = eam_potentials

  @property
  def type(self):
    return "EAM"

  @property
  def nrho(self):
    return self._nrho

  @property
  def cutoff_rho(self):
    return self._cutoff_rho

  @property
  def drho(self):
    return (self.cutoff_rho / float(self.nrho-1))

  @property
  def eam_potentials(self):
    return self._eam_potentials


class SetFL_EAMTabulation(_EAMTabulationAbstractbase):
  """Class for tabulating setfl formatted embedded atom potentials suitable
  for use with LAMMPS' pair_style eam/alloy"""

  def __init__(self, potentials, eam_potentials, cutoff, nr, cutoff_rho, nrho):
    """Instantiate class for tabulation of setfl formatted embedded atom potential tables.

    :params potentials: List of atsim.potentials.Potential objects.
    :params eam_potentials: List of `atsim.potentials.EAMPotential` instances.
    :params cutoff: Maximum separation to be tabulated.
    :params nr: Number of points to be used in tabulation
    :params cutoff_rho: Density cutoff.
    :params nrho: Number of points to be used when discretising density range during EAM tabulation"""
    super(SetFL_EAMTabulation, self).__init__(potentials, eam_potentials, cutoff, nr, cutoff_rho, nrho, "setfl")

  def write(self, fp):
    """Write the tabulation to the file object `fp`.

    :param fp: File object into which data should be written."""
    writeSetFL(
      self.nrho, self.drho, 
      self.nr, self.dr,
      self.eam_potentials,
      self.potentials,
      out = fp)

class SetFL_FS_EAMTabulation(_EAMTabulationAbstractbase):
  """Class for tabulating setfl Finnis-Sinclair formatted embedded atom potentials suitable
  for use with LAMMPS' pair_style eam/fs"""

  def __init__(self, potentials, eam_potentials, cutoff, nr, cutoff_rho, nrho):
    """Instantiate class for tabulation of setfl Finnis-Sinclair formatted embedded atom potential tables.

    :params potentials: List of atsim.potentials.Potential objects.
    :params eam_potentials: List of `atsim.potentials.EAMPotential` instances.
    :params cutoff: Maximum separation to be tabulated.
    :params nr: Number of points to be used in tabulation
    :params cutoff_rho: Density cutoff.
    :params nrho: Number of points to be used when discretising density range during EAM tabulation"""
    super(SetFL_FS_EAMTabulation, self).__init__(potentials, eam_potentials, cutoff, nr, cutoff_rho, nrho, "setfl_fs")

  def write(self, fp):
    """Write the tabulation to the file object `fp`.

    :param fp: File object into which data should be written."""
    writeSetFLFinnisSinclair(
      self.nrho, self.drho, 
      self.nr, self.dr,
      self.eam_potentials,
      self.potentials,
      out = fp)

class TABEAM_EAMTabulation(_EAMTabulationAbstractbase):
  """Class for tabulating TABEAM formatted embedded atom potentials for the DL_POLY code."""

  def __init__(self, potentials, eam_potentials, cutoff, nr, cutoff_rho, nrho):
    """Instantiate class for tabulation of DL_POLY TABEAM formatted embedded atom potential tables.

    :params potentials: List of atsim.potentials.Potential objects.
    :params eam_potentials: List of `atsim.potentials.EAMPotential` instances.
    :params cutoff: Maximum separation to be tabulated.
    :params nr: Number of points to be used in tabulation
    :params cutoff_rho: Density cutoff.
    :params nrho: Number of points to be used when discretising density range during EAM tabulation"""
    super(TABEAM_EAMTabulation, self).__init__(potentials, eam_potentials, cutoff, nr, cutoff_rho, nrho, "DL_POLY_EAM")

  def write(self, fp):
    """Write the tabulation to the file object `fp`.

    :param fp: File object into which data should be written."""
    writeTABEAM(
      self.nrho, self.drho, 
      self.nr, self.dr,
      self.eam_potentials,
      self.potentials,
      out = fp)

class TABEAM_FinnisSinclair_EAMTabulation(_EAMTabulationAbstractbase):
  """Class for tabulating EEAM TABEAM formatted Finnis-Sinclair style embedded atom potentials for the DL_POLY code."""

  def __init__(self, potentials, eam_potentials, cutoff, nr, cutoff_rho, nrho):
    """Instantiate class for tabulation of DL_POLY EEAM TABEAM formatted Finnis-Sinclair embedded atom potential tables.

    :params potentials: List of atsim.potentials.Potential objects.
    :params eam_potentials: List of `atsim.potentials.EAMPotential` instances.
    :params cutoff: Maximum separation to be tabulated.
    :params nr: Number of points to be used in tabulation
    :params cutoff_rho: Density cutoff.
    :params nrho: Number of points to be used when discretising density range during EAM tabulation"""
    super(TABEAM_FinnisSinclair_EAMTabulation, self).__init__(potentials, eam_potentials, cutoff, nr, cutoff_rho, nrho, "DL_POLY_EAM_fs")

  def write(self, fp):
    """Write the tabulation to the file object `fp`.

    :param fp: File object into which data should be written."""
    writeTABEAMFinnisSinclair(
      self.nrho, self.drho, 
      self.nr, self.dr,
      self.eam_potentials,
      self.potentials,
      out = fp)

