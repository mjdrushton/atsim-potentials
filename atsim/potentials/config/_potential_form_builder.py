import logging

from .._multi_range_potential_form import Multi_Range_Potential_Form

class Potential_Form_Builder(object):
  """Class that instantiates potential function callables from PotentialFormInstanceTuple"""

  def __init__(self, potential_form_registry, modifier_registry):
    """:param potential_form_registry: Instance of `atsim.potentials.config._potential_form_registry.Potential_Form_Registry`
       :param modifier_registry: Instance of `atsim.potentials.config._modifier_registry.Modifier_Registry`"""
    self.potential_form_registry = potential_form_registry
    self.modifier_registry = modifier_registry
    
  def _make_multi_range_tuple(self, pform_instance):
    """Convert PotentialFormInstanceTuple into Multi_Range_Potential_Form.Range_Defn_Tuple"""

    # Is pform_instance a modifier or a potential_instance?
    if hasattr(pform_instance, "modifier"):
      pform_factory = self.modifier_registry[pform_instance.modifier]
      pform  = pform_factory(pform_instance.potential_forms, self)
    else:
      pform_factory = self.potential_form_registry[pform_instance.potential_form]
      params = pform_instance.parameters
      pform = pform_factory(*params)

    if pform_instance.start:
      start = pform_instance.start.start
      range_type = pform_instance.start.range_type
    else:
      start = float("-inf")
      range_type = ">="

    mr_tuple = Multi_Range_Potential_Form.Range_Defn_Tuple(range_type, start, pform)
    return mr_tuple

  def create_potential_function(self, potential_form_instance):
    """Create a callable for a given potential form tuple.

       :param potential_form_instance: PotentialFormInstanceTuple (or PotentialModifierTuple) defining potential form callable.
       :return: Single parameter potential function"""   
    logger = logging.getLogger(__name__).getChild("Potential_Form_Builder.create_potential")
    logger.debug("Creating potential object for potential form instance: {}".format(potential_form_instance))
    
    # Parametrise the potential form to create a potential function
    tuples = [self._make_multi_range_tuple(potential_form_instance)]
    n = potential_form_instance.next
    while n:
      tuples.append(self._make_multi_range_tuple(n))
      n = n.next
    pot_func = Multi_Range_Potential_Form(*tuples)
    return pot_func