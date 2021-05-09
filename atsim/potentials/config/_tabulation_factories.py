import logging

import collections

try:
  import collections.abc
  Mapping = collections.abc.Mapping
except ImportError:
  Mapping = collections.Mapping

from ..pair_tabulation import DLPoly_PairTabulation, LAMMPS_PairTabulation, GULP_PairTabulation, Excel_PairTabulation
from ..eam_tabulation  import SetFL_EAMTabulation, SetFL_FS_EAMTabulation, TABEAM_EAMTabulation, TABEAM_FinnisSinclair_EAMTabulation, Excel_EAMTabulation, Excel_FinnisSinclair_EAMTabulation, ADP_EAMTabulation

from ._common import ConfigurationException
from ._potential_form_registry import Potential_Form_Registry
from ._modifier_registry import Modifier_Registry
from ._pair_potential_builder import Pair_Potential_Builder, Pair_Potentials_From_Tuples_Builder
from ._eam_potential_builder import EAM_Potential_Builder, EAM_Potential_Builder_FS
from ..referencedata import Reference_Data

RCutoffTuple = collections.namedtuple('RCutoffTuple', ['cutoff', 'nr'])
R_Rho_CutoffTuple = collections.namedtuple('R_Rho_CutoffTuple', ['cutoff', 'nr', 'cutoff_rho', 'nrho'])

def _create_pair_objects(potential_form_registry, modifier_registry, cp):
  potbuilder = Pair_Potential_Builder(cp, potential_form_registry, modifier_registry)
  pots = potbuilder.potentials
  return pots

class PairTabulationFactory(object):
  """Factory class for Tabulation objects such as those in the atsim.potentials.pair_tabulation module.

  The role of this class is to extract data from an atsim.potentials.config.ConfigParser and together
  with other builder classes, create atsim.potentials.Potential objects and extract the cutoff values that
  can be used to instantiate a Tabulation object.
  
  The key method for clients of this object is .create_tabulation().

  Other public methods are provided to allow customisation in sub-classes.
  
  """

  def __init__(self, tabulation_target, tabulation_class):
    self.tabulation_target = tabulation_target
    self.tabulation_class = tabulation_class
    self.tabulation_type = "pair-potential"

  def _log_cutoffs(self, logger, r_cutoff, **kwargs):
    logger.info("  * cutoff: {}".format(r_cutoff.cutoff))
    logger.info("  * nr: {}".format(r_cutoff.nr))

  def _log_extra_tabulation_details(self, logger, r_cutoff, potobjs, **kwargs):
    pass

  def _log_pair_potentials(self, logger, potobjs, **kwargs):
    logger.info("  * Tabulation will contain the following pair-potentials:")
    for pot in potobjs:
      logger.info("      + {}-{}".format(pot.speciesA, pot.speciesB))

  def _log_more(self, logger, r_cutoff, potobjs, **kwargs):
    pass

  def _log_tabulation_details(self, r_cutoff, potobjs, **kwargs):
    logger = logging.getLogger(__name__).getChild("PairTabulationFactory._log_tabulation_details")

    logger.info("Creating {tabulation_type} tabulation with following properties:".format(tabulation_type = self.tabulation_type))
    logger.info("  * tabulation-target : {}".format(self.tabulation_target))
    self._log_cutoffs(logger, r_cutoff, **kwargs)
    self._log_extra_tabulation_details(logger,r_cutoff, potobjs, **kwargs)
    self._log_pair_potentials(logger, potobjs, **kwargs)
    self._log_more(logger, r_cutoff, potobjs, **kwargs)

  def extract_cutoffs(self, cp):
    """Called by .create_tabulation() to determine the maximum separation cutoff for the tabulation.

    :param cp: atsim.potentials.config.ConfigParser instance representing parsed configuration file.
    :return: RCutoffTuple (or object providing .cutoff and .nr properties)"""

    logger = logging.getLogger(__name__).getChild("PairTabulationFactory.extract_cutoffs")

    # Get cutoff and gridpoints
    if cp.tabulation.cutoff is None:
      cutoff = 10.0
      logger.warning("Cutoff not specified using a default of {}".format(cutoff))
    else:
      cutoff = cp.tabulation.cutoff
    
    if cp.tabulation.nr is None:
      nr = 1001
      logger.warning("nr not specified using a default of {}".format(nr))
    else:
      nr = cp.tabulation.nr

    return RCutoffTuple(cutoff, nr)

  def extract_potential_objects(self, cp, potential_form_registry, modifier_registry):
    """Method called by .create_tabulation() to create a list of atsim.potentials.Potential objects.
    
    :param cp: atsim.potentials.config.ConfigParser instance representing parsed configuration file.
    :param potential_form_registry: atsim.potentials.config.Potential_Form_Registry object.
    :param modifier_registry: atsim.potentials.config.Modifier_Registry instance.

    :return: List of atsim.potentials.Potential objects.
    """

    potobjs = _create_pair_objects(potential_form_registry,modifier_registry, cp)
    return potobjs

  def extract_tabulation_args(self, cp, r_cutoff, potobjs, potential_form_registry, modifier_registry):
    """Extract and prepare arguments to be passed to Tabulation class constructor.

    Returns a list of arguments suitable for passing to Tabulation class constructor. Typically these are: 
    ``[potobjs, cutoff, nr]`` 

    Where:
      * ``potobjs`` is a list of :class:`atsim.potentials.Potential` instances.
      * ``cutoff`` is a float giving the maximum separation for the tabulation.
      * ``nr`` the number of rows to be used in the tabulation.

    :param cp: atsim.potentials.config.ConfigParser instance representing parsed configuration file.
    :param r_cutoff: Float giving the maximum separation for tabulation (as extracted from configuration).
    :param potobjs: List of atsim.potentials.Potential objects to be tabulated.
    :param potential_form_registry: atsim.potentials.config.Potential_Form_Registry object.
    :param modifier_registry: atsim.potentials.config.Modifier_Registry instance.
    
    :return: List containing Tabulation class constructor arguments. """

    return [potobjs, r_cutoff.cutoff, r_cutoff.nr]

  def create_tabulation(self,cp):
    # logger = logging.getLogger(__name__).getChild("PairTabulationFactory.create_tabulation")
    r_cutoff = self.extract_cutoffs(cp)

    # Get pair potentials
    potential_form_registry = Potential_Form_Registry(cp, register_standard = True, register_pymath_functions = True)
    modifier_registry = Modifier_Registry()
    
    potobjs = self.extract_potential_objects(cp, potential_form_registry, modifier_registry)
    self._log_tabulation_details(r_cutoff, potobjs)

    tabulationargs = self.extract_tabulation_args(cp, r_cutoff, potobjs, potential_form_registry, modifier_registry)
    tabulation = self.tabulation_class(*tabulationargs)
    return tabulation


class EAMTabulationFactory(PairTabulationFactory):
  """Tabulation factory for setfl (LAMMPS eam/alloy) potentials"""

  def __init__(self, tabulation_target, tabulation_class, eam_builder_class = EAM_Potential_Builder):
    super(EAMTabulationFactory, self).__init__(tabulation_target, tabulation_class)
    self.eam_builder_class = eam_builder_class
    self.tabulation_type = "EAM"

  def extract_cutoffs(self, cp):
    logger = logging.getLogger(__name__).getChild("EAMTabulationFactory._get_cutoffs")
    r_cutoff = super(EAMTabulationFactory, self).extract_cutoffs(cp)

    # Get cutoff and gridpoints
    if cp.tabulation.cutoff_rho is None:
      cutoff_rho = 100.0
      logger.warning("cutoff_rho not specified using a default of {}".format(cutoff_rho))
    else:
      cutoff_rho = cp.tabulation.cutoff_rho
    
    if cp.tabulation.nrho is None:
      nrho = 1001
      logger.warning("nrho not specified using a default of {}".format(nrho))
    else:
      nrho = cp.tabulation.nrho
    return R_Rho_CutoffTuple(r_cutoff.cutoff, r_cutoff.nr, cutoff_rho, nrho)

  def _create_reference_data(self, cp):
    extra_rd = cp.species
    rd = Reference_Data(extra_rd)
    return rd

  def extract_tabulation_args(self, cp, r_cutoff, potobjs, potential_form_registry, modifier_registry):
    logger = logging.getLogger(__name__).getChild("EAMTabulationFactory._make_tabulation_args")
    rd = self._create_reference_data(cp)
    eam_builder =  self.eam_builder_class(cp, potential_form_registry, modifier_registry, rd)
    eam_potentials = eam_builder.eam_potentials

    logger.info("  * Tabulation will contain following EAM species:")
    any_dict = False
    for eam_potential in eam_potentials:
      logger.info("      + {}".format(eam_potential.species))
      if isinstance(eam_potential.electronDensityFunction, Mapping):
        any_dict = True
    
    if any_dict:
      logger.info("  * Tabulation will contain the folowing density functions:")
      for eam_potential in eam_potentials:
        if isinstance(eam_potential.electronDensityFunction, Mapping):
          other_species = sorted(eam_potential.electronDensityFunction.keys())
          for s in other_species:
            logger.info("      + {}->{}".format(eam_potential.species, s))
        else:
            logger.info("      + {}".format(eam_potential.species))

    args = [potobjs, eam_potentials, 
            r_cutoff.cutoff, r_cutoff.nr, 
            r_cutoff.cutoff_rho, r_cutoff.nrho]
    return args

  def _log_cutoffs(self, logger, r_cutoff, **kwargs):
    super(EAMTabulationFactory, self)._log_cutoffs(logger, r_cutoff, **kwargs)
    logger.info("  * cutoff_rho: {}".format(r_cutoff.cutoff_rho))
    logger.info("  * nrho: {}".format(r_cutoff.nrho))

class DLPOLY_PairTabulationFactory(PairTabulationFactory):
  """PairTabulationFactory that performs extra checks required by DL_POLY TABLEs and raises ConfigException if errors are found"""

  def extract_cutoffs(self, cp):
    cutoffs = super(DLPOLY_PairTabulationFactory, self).extract_cutoffs(cp)
    if cutoffs.nr % 4 != 0:
      raise ConfigurationException("The number of rows in a DL_POLY TABLE file needs to be divisible by 4. Number of rows specified = {} ".format(cutoffs.nr))
    return cutoffs

class ADP_EAMTabulationFactory(EAMTabulationFactory):
  """EAMTabulationFactory which creates the additional dipole and quadrupole objects 
  required by the ADP EAM extension"""

  def _extract_pots(self, cp, pfr, mr, section_name):
    tuples = cp.parse_pair_like(section_name)
    pot_builder = Pair_Potentials_From_Tuples_Builder(tuples, pfr, mr, section_name )
    return pot_builder.potentials

  def extract_dipoles(self, cp, potential_form_registry, modifier_registry):
    return self._extract_pots(cp, potential_form_registry, modifier_registry, "EAM-ADP-Dipole")

  def extract_quadrupoles(self, cp, potential_form_registry, modifier_registry):
    return self._extract_pots(cp, potential_form_registry, modifier_registry, "EAM-ADP-Quadrupole")

  def extract_tabulation_args(self, cp, r_cutoff, potobjs, potential_form_registry, modifier_registry):
    (potobjs, eam_potentials, 
    cutoff_r, nr, 
    cutoff_rho, nrho) = super().extract_tabulation_args(cp, r_cutoff, potobjs, potential_form_registry, modifier_registry)

    dipoles = self.extract_dipoles(cp, potential_form_registry, modifier_registry)
    quadrupoles = self.extract_quadrupoles(cp, potential_form_registry, modifier_registry)
    args = [potobjs, eam_potentials, dipoles, quadrupoles, cutoff_r, nr, cutoff_rho, nrho]
    return args

"""Target name to factory objects"""
TABULATION_FACTORIES = {
  "LAMMPS"       :  PairTabulationFactory("LAMMPS", LAMMPS_PairTabulation),
  "DLPOLY"       :  DLPOLY_PairTabulationFactory("DLPOLY", DLPoly_PairTabulation),
  "GULP"         :  PairTabulationFactory("GULP", GULP_PairTabulation),
  "excel"        :  PairTabulationFactory("excel", Excel_PairTabulation),
  "setfl"        :  EAMTabulationFactory("setfl/lammps_eam_alloy", SetFL_EAMTabulation),    
  "setfl_fs"     :  EAMTabulationFactory("setfl/lammps_eam_fs", SetFL_FS_EAMTabulation, EAM_Potential_Builder_FS),
  "DL_POLY_EAM"  :  EAMTabulationFactory("DL_POLY_EAM", TABEAM_EAMTabulation),               
  "DL_POLY_EAM_fs" :  EAMTabulationFactory("DL_POLY_EAM_fs", TABEAM_FinnisSinclair_EAMTabulation, EAM_Potential_Builder_FS),
  "excel_eam"    :  EAMTabulationFactory("excel_eam", Excel_EAMTabulation),
  "excel_eam_fs"    :  EAMTabulationFactory("excel_eam_fs", Excel_FinnisSinclair_EAMTabulation, EAM_Potential_Builder_FS),
  "eam_adp" : ADP_EAMTabulationFactory("eam_adp", ADP_EAMTabulation)
}
