import logging

import collections

from .._pair_tabulation import DLPoly_PairTabulation, LAMMPS_PairTabulation
from .._eam_tabulation  import SetFL_EAMTabulation, SetFL_FS_EAMTabulation, TABEAM_EAMTabulation, TABEAM_FinnisSinclair_EAMTabulation

from ._potential_form_registry import Potential_Form_Registry
from ._modifier_registry import Modifier_Registry
from ._pair_potential_builder import Pair_Potential_Builder
from ._eam_potential_builder import EAM_Potential_Builder, EAM_Potential_Builder_FS
from ..referencedata import Reference_Data

RCutoffTuple = collections.namedtuple('RCutoffTuple', ['cutoff', 'nr'])
R_Rho_CutoffTuple = collections.namedtuple('R_Rho_CutoffTuple', ['cutoff', 'nr', 'cutoff_rho', 'nrho'])

def _create_pair_objects(potential_form_registry, modifier_registry, cp):
  potbuilder = Pair_Potential_Builder(cp, potential_form_registry, modifier_registry)
  pots = potbuilder.potentials
  return pots

class PairTabulationFactory(object):

  def __init__(self, tabulation_target, tabulation_class):
    self.tabulation_target = tabulation_target
    self.tabulation_class = tabulation_class
    self.tabulation_type = "pair-potential"

  def _get_cutoffs(self, cp):
    logger = logging.getLogger(__name__).getChild("PairTabulationFactory._get_cutoffs")

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

  def _make_tabulation_args(self, cp, r_cutoff, potobjs, potential_form_registry, modifier_registry):
    return [potobjs, r_cutoff.cutoff, r_cutoff.nr]

  def create_tabulation(self,cp):
    logger = logging.getLogger(__name__).getChild("PairTabulationFactory.create_tabulation")
    r_cutoff = self._get_cutoffs(cp)

    # Get pair potentials
    potential_form_registry = Potential_Form_Registry(cp, True)
    modifier_registry = Modifier_Registry()
    potobjs = _create_pair_objects(potential_form_registry,modifier_registry, cp)
    self._log_tabulation_details(r_cutoff, potobjs)

    tabulationargs = self._make_tabulation_args(cp, r_cutoff, potobjs, potential_form_registry, modifier_registry)

    tabulation = self.tabulation_class(*tabulationargs)
    return tabulation


class EAMTabulationFactory(PairTabulationFactory):
  """Tabulation factory for setfl (LAMMPS eam/alloy) potentials"""

  def __init__(self, tabulation_target, tabulation_class, eam_builder_class = EAM_Potential_Builder):
    super(EAMTabulationFactory, self).__init__(tabulation_target, tabulation_class)
    self.eam_builder_class = eam_builder_class
    self.tabulation_type = "EAM"

  def _get_cutoffs(self, cp):
    logger = logging.getLogger(__name__).getChild("EAMTabulationFactory._get_cutoffs")
    r_cutoff = super(EAMTabulationFactory, self)._get_cutoffs(cp)

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

  def _make_tabulation_args(self, cp, r_cutoff, potobjs, potential_form_registry, modifier_registry):
    logger = logging.getLogger(__name__).getChild("EAMTabulationFactory._make_tabulation_args")
    rd = self._create_reference_data(cp)
    eam_builder =  self.eam_builder_class(cp, potential_form_registry, modifier_registry, rd)
    eam_potentials = eam_builder.eam_potentials

    logger.info("  * Tabulation will contain following EAM species:")
    any_dict = False
    for eam_potential in eam_potentials:
      logger.info("      + {}".format(eam_potential.species))
      if isinstance(eam_potential.electronDensityFunction, collections.Mapping):
        any_dict = True
    
    if any_dict:
      logger.info("  * Tabulation will contain the folowing density functions:")
      for eam_potential in eam_potentials:
        if isinstance(eam_potential.electronDensityFunction, collections.Mapping):
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

"""Target name to factory objects"""
TABULATION_FACTORIES = {
  "LAMMPS"       :  PairTabulationFactory("LAMMPS", LAMMPS_PairTabulation),                 
  "DLPOLY"       :  PairTabulationFactory("DLPOLY", DLPoly_PairTabulation),                 
  "setfl"        :  EAMTabulationFactory("setfl/lammps_eam_alloy", SetFL_EAMTabulation),    
  "setfl_fs"     :  EAMTabulationFactory("setfl/lammps_eam_fs", SetFL_FS_EAMTabulation, EAM_Potential_Builder_FS),
  "DL_POLY_EAM"  :  EAMTabulationFactory("DL_POLY_EAM", TABEAM_EAMTabulation),               
  "DL_POLY_EAM_fs" :  EAMTabulationFactory("DL_POLY_EAM_fs", TABEAM_FinnisSinclair_EAMTabulation, EAM_Potential_Builder_FS)               
}
