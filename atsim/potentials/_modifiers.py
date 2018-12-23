import functools
import logging

from .config._common import ConfigurationException, MultiRangeDefinitionTuple

from atsim.potentials import plus

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

  sum_mod = functools.reduce(plus, pot_callables)

  return sum_mod

@modifier
def spline(potential_forms, potential_form_builder):
  """Modifier that smoothly splines between two potential forms by linking them with an intermediate spline.

  The `potential_forms` list must contain a single `PotentialFormInstanceTuple` entre. The tuple
  must only define three sections - potential form A -> exp_spline -> potential form B. The spline is defined for the region where exp_spline starts and potential form B starts. 

  As a configuration string this might be define as:

  `>0 as.zbl 14 8 >=0.8 exp_spline >=1.4 as.buck 180003 0.3 32.0`

  Which would create a `zbl` and Buckingham potential connected by a spline when `r` is between 0.8 and 1.4.

  At present the potential form label for the spline must be `exp_spline`.    

  :param potential_forms: List of tuples that can be passed to `atsim.potentials.config._potential_form_builder.Potential_Form_Builder.create_potential_function()` to create potential callables.
  :param potential_form_builder: `atsim.potentials.config._potential_form_builder.Potential_Form_Builder` used to create potential instances.

  :returns: Potential callable representing splined potentials."""

  if len(potential_forms) != 1:
    raise ConfigurationException("spline modifier only accepts a single multi range potential definition as its argument")

  logger = logging.getLogger(__name__).getChild("spline")

  # Extract the ranges and potential forms.
  pform = potential_forms[0]
  if pform.next is None:
    raise ConfigurationException("spline modifier requires three sub-potentials to be defined only one specified.")

  # ... set the 'next' attribute to None
  tdict = pform._asdict()
  tdict['next'] = None
  pot1 = pform.__class__(**tdict)

  tdict = pform.next._asdict()
  tdict['next'] = None
  pot2 = pform.next.__class__(**tdict)

  if pot2.potential_form != "exp_spline":
    raise ConfigurationException("spline modifier only accepts 'exp_spline' middle potential form. '{}' was found instead".format(pot2.potential_form))

  if pot2.parameters:
    raise ConfigurationException("spline modifier 'exp_spline' middle potential form does not take any parameters. The following parameters were specified: {}".format(pot2.parameters))

  if pform.next.next is None:
    raise ConfigurationException("spline modifier requires three sub-potentials to be defined only two specified.")

  tdict = pform.next.next._asdict()
  tdict['next'] = None
  # ... we need to be able to evaluate the potential at the start point
  #     this allows values to be calculated at any point (as long as the potential's well behaved there)
  pot3_start =tdict['start'].start
  tdict['start'] = MultiRangeDefinitionTuple(u">", float("-inf"))
  pot3 = pform.next.next.__class__(**tdict)

  if not pform.next.next.next is None:
    raise ConfigurationException("spline modifier requires three sub-potentials to be defined more than three have been given.")

  from ._spline import SplinePotential

  # Determine detachment point
  if not pot1.start.start < pot2.start.start:
    raise ConfigurationException("spline modifier range error. Start of 1st potential should be less than start of 2nd ! {} < {}".format(pot1.start.start, pot2.start.start))

  if not pot2.start.start < pot3_start:
    raise ConfigurationException("spline modifier range error. Start of 2nd potential (exp_spline) should be less than start of 3rd ! {} < {}".format(pot2.start.start, pot3_start))

  detach_point = pform.next.start.start
  attach_point = pform.next.next.start.start

  pot1_func = potential_form_builder.create_potential_function(pot1)
  pot3_func = potential_form_builder.create_potential_function(pot3)

  logger.debug("spline modifier: connecting '{}' with exp_spline to '{}' in range {} to {}".format(pot1.potential_form, pot2.potential_form, detach_point, attach_point))

  # Now build the spline object
  spot_obj = SplinePotential(pot1_func, pot3_func, detach_point, attach_point)
  return spot_obj


