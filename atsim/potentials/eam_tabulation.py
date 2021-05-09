
from .pair_tabulation import PairTabulation_AbstractBase, Excel_PairTabulation, _r_value_iterator

from ._lammpsWriteEAM import writeSetFL, writeSetFLFinnisSinclair, _writeSetFLPairPots
from ._dlpoly_writeTABEAM import writeTABEAM, writeTABEAMFinnisSinclair

def _rho_value_iterator(tabulation):
  #for n in range(tabulation.nr+1):
  for n in range(tabulation.nrho):
    yield float(n)* tabulation.cutoff_rho / (float(tabulation.nrho) -1)


class _EAMTabulationAbstractbase(PairTabulation_AbstractBase):
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


class Excel_EAMTabulation(_EAMTabulationAbstractbase):
  """Class for dumping EAM model into a spreadsheet"""

  _excel_tab_name = "excel_eam"

  def __init__(self, potentials, eam_potentials, cutoff, nr, cutoff_rho, nrho):
    """Instantiate class for tabulation of setfl formatted embedded atom potential tables.

    :params potentials: List of atsim.potentials.Potential objects.
    :params eam_potentials: List of `atsim.potentials.EAMPotential` instances.
    :params cutoff: Maximum separation to be tabulated.
    :params nr: Number of points to be used in tabulation
    :params cutoff_rho: Density cutoff.
    :params nrho: Number of points to be used when discretising density range during EAM tabulation"""
    super(Excel_EAMTabulation, self).__init__(potentials, eam_potentials, cutoff, nr, cutoff_rho, nrho, self._excel_tab_name )
    self._inner_tabulation = None 

  def _build_workbook(self):
    self._inner_tabulation = Excel_PairTabulation(self.potentials, self.cutoff, self.nr)
    wb = self._inner_tabulation.workbook
    self._add_sheets(wb)

  def _add_sheets(self, wb):
    self._add_eam_density(wb)
    self._add_eam_embed(wb)


  def _add_eam_density(self, wb):
    ws = wb.create_sheet("EAM-Density")

    # Get sorted list of potentials
    pot_dict = {}
    for p in self.eam_potentials:
      k = p.species
      v = p.electronDensityFunction
      pot_dict[k] = v
    column_heads = sorted(pot_dict.keys())
    self._inner_tabulation._populate_worksheet(ws, "r", _r_value_iterator(self), column_heads, pot_dict )

  def _add_eam_embed(self, wb):
    ws = wb.create_sheet("EAM-Embed")

    # Get sorted list of potentials
    pot_dict = {}
    for p in self.eam_potentials:
      k = p.species
      v = p.embeddingFunction
      pot_dict[k] = v
    column_heads = sorted(pot_dict.keys())
    self._inner_tabulation._populate_worksheet(ws, "rho", _rho_value_iterator(self), column_heads, pot_dict )

  @property
  def workbook(self):
    if self._inner_tabulation is None:
      self._build_workbook()
    return self._inner_tabulation.workbook

  def write(self, fp):
    """Write the tabulation to the file object `fp`.

    :param fp: File object into which data should be written."""
    wb = self.workbook
    self._inner_tabulation.write(fp)

  @classmethod
  def open_fp(cls, filename):
    """Creates a file object with a given path suitable for writing potential data to.

    :param filename: Filename of output file object.

    :return: File object suitable for passing to write() method"""
    return Excel_PairTabulation.open_fp(filename)


class Excel_FinnisSinclair_EAMTabulation(Excel_EAMTabulation):
  """Class for dumping EAM model into a spreadsheet"""

  _excel_tab_name = "excel_eam_fs"
    

  def _add_eam_density(self, wb):
    ws = wb.create_sheet("EAM-Density")

    # Get sorted list of potentials
    pot_dict = {}
    for p in self.eam_potentials:
      species_f = p.species
      v = p.electronDensityFunction
      for species_t, func in p.electronDensityFunction.items():
        k = "{}->{}".format(species_f, species_t)
        pot_dict[k] = func

    column_heads = sorted(pot_dict.keys())
    self._inner_tabulation._populate_worksheet(ws, "r", _r_value_iterator(self), column_heads, pot_dict )


class ADP_EAMTabulation(SetFL_EAMTabulation):
  """Class for tabulating setfl formatted embedded atom potentials with the ADP, angular dependent extension,
  suitable for use with LAMMPS' pair_style adp"""

  def __init__(self, potentials, eam_potentials, dipole_potentials, quadrupole_potentials, cutoff, nr, cutoff_rho, nrho):
    """Instantiate class for tabulation of setfl formatted embedded atom potential tables.

    :params potentials: List of atsim.potentials.Potential objects.
    :params eam_potentials: List of `atsim.potentials.EAMPotential` instances.
    :params dipole_potentials: List atsim.potentials.Potential objects giving ADP dipole functions.
    :params quadrupole_potentials: List atsim.potentials.Potential objects giving ADP quadrupole functions.
    :params cutoff: Maximum separation to be tabulated.
    :params nr: Number of points to be used in tabulation
    :params cutoff_rho: Density cutoff.
    :params nrho: Number of points to be used when discretising density range during EAM tabulation"""
    super(SetFL_EAMTabulation, self).__init__(potentials, eam_potentials, cutoff, nr, cutoff_rho, nrho, "eam_adp")
    self.dipole_potentials = dipole_potentials
    self.quadrupole_potentials = quadrupole_potentials

  def write(self, fp):
    """Write the tabulation to the file object `fp`.

    :param fp: File object into which data should be written."""
    writeSetFL(
      self.nrho, self.drho, 
      self.nr, self.dr,
      self.eam_potentials,
      self.potentials,
      out = fp)

    self._write_dipole(fp)
    self._write_quadrupole(fp)


  def _write_dipole(self, fp):
    _writeSetFLPairPots(self.nr, self.dr, self.eam_potentials, self.dipole_potentials, fp, scale_r=False)

  def _write_quadrupole(self, fp):
    _writeSetFLPairPots(self.nr, self.dr, self.eam_potentials, self.quadrupole_potentials, fp, scale_r=False)