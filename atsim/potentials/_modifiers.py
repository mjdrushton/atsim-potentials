import functools
import logging

def modifier(func):
  func.is_modifier = True
  return func

def is_modifier(obj):
  return hasattr(obj, "is_modifier") and obj.is_modifier

# As we're about to shadow it - save the builtin sum to _sum so we can still access it.
_sum = sum

@modifier
def sum(potential_forms, potential_form_builder):
  """Modifier that sums all the potential instances given as arguments.

  :param potential_forms: List of tuples that can be passed to `atsim.potentials.config._potential_form_builder.Potential_Form_Builder.create_potential_function()` to create potential callables.
  :param potential_form_builder: `atsim.potentials.config._potential_form_builder.Potential_Form_Builder` used to create potential instances.

  :returns: Potential callable sums that values from a number of potential instances."""
  logger = logging.getLogger(__name__).getChild("add")

  logger.debug("Creating 'add' modifier for:")
  pot_callables = []
  for i,pfi in enumerate(potential_forms):
    logger.debug("  {}: {}".format(i+1, pfi))
    pot_callable = potential_form_builder.create_potential_function(pfi)
    pot_callables.append(pot_callable)

  def sum_func(r):
    values = [pc(r) for pc in pot_callables]
    summed = _sum(values)
    return summed
  return sum_func