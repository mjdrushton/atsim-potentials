import logging

from .._multi_range_potential_form import Multi_Range_Potential_Form

class Potential_Form_Builder(object):
  """Class that instantiates potential function callables from PotentialFormInstanceTuple"""


  def __init__(self, potential_form_registry):
    """:param potential_form_registry: Instance of `atsim.potentials.config._potential_form_registry.Potential_Form_Registry`"""
    self.potential_form_registry = potential_form_registry

  def _make_multi_range_tuple(self, pform_instance):
    """Convert PotentialFormInstanceTuple into Multi_Range_Potential_Form.Range_Defn_Tuple"""
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

       :param potential_form_instance: PotentialFormInstanceTuple defining potential form callable.
       :return: Single parameter potential function"""   
    logger = logging.getLogger(__name__).getChild("Potential_Form_Builder.create_potential")
    logger.debug("Creating potential object for potential form instance: {}".format(potential_form_instance))
    potential_form = self.potential_form_registry[potential_form_instance.potential_form]
    
    # Parametrise the potential form to create a potential function
    tuples = [self._make_multi_range_tuple(potential_form_instance)]
    n = potential_form_instance.next
    while n:
      tuples.append(self._make_multi_range_tuple(n))
      n = n.next
    pot_func = Multi_Range_Potential_Form(*tuples)
    return pot_func
    

# class Multi_Range_Potential_Form_Builder(object):
#   """Class that instantiates composite potential functions callables.
#   These are composed of multiple sub-functions that each describe a particular range of input values"""

#   def __init__(self, potential_form_builder):
#     """:param potential_form_builder: Instance of Potential_Form_Builder used to create sub functions"""
#     self.potential_form_builder = potential_form_builder

#   def create_potential_function(self, potential_form_instance):
#     """Create a callable for a given potential form tuple.

#        :param potential_form_instance: PotentialFormInstanceTuple defining potential form callable.
#        :return: Single parameter potential function"""   
#     logger = logging.getLogger(__name__).getChild("Multi_Range_Potential_Form_Builder.create_potential")
#     logger.debug("Creating potential object for potential form instance: {}".format(potential_form_instance))



