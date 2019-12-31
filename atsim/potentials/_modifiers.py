import functools
import logging

from .config._common import ConfigurationException, MultiRangeDefinitionTuple

from atsim.potentials import plus

from .spline import Custom_SplinePotential, Spline_Point, Exp_Spline, Buck4_Spline


def modifier(func):
  func.is_modifier = True
  return func

def is_modifier(obj):
  return hasattr(obj, "is_modifier") and obj.is_modifier

# As we're about to shadow it - save the builtin sum to _sum so we can still access it.
_sum = sum

def _modifier_from_func_reduce(logger_name, func, potential_forms, potential_form_builder):
  logger = logging.getLogger(__name__).getChild(logger_name)

  logger.debug("Creating '{}' modifier for:".format(logger_name))
  pot_callables = []
  for i,pfi in enumerate(potential_forms):
    logger.debug("  {}: {}".format(i+1, pfi))
    pot_callable = potential_form_builder.create_potential_function(pfi)
    pot_callables.append(pot_callable)

  mod = functools.reduce(func, pot_callables)

  return mod


@modifier
def sum(potential_forms, potential_form_builder):
  """Modifier that sums all the potential instances given as arguments.

  :param potential_forms: List of tuples that can be passed to `atsim.potentials.config._potential_form_builder.Potential_Form_Builder.create_potential_function()` to create potential callables.
  :param potential_form_builder: `atsim.potentials.config._potential_form_builder.Potential_Form_Builder` used to create potential instances.

  :returns: Potential callable sums that values from a number of potential instances."""
  mod = _modifier_from_func_reduce("sum", plus, potential_forms, potential_form_builder)
  return mod

@modifier
def product(potential_forms, potential_form_builder):
  """Modifier that takes the product of the potential instances given as arguments.

  :param potential_forms: List of tuples that can be passed to `atsim.potentials.config._potential_form_builder.Potential_Form_Builder.create_potential_function()` to create potential callables.
  :param potential_form_builder: `atsim.potentials.config._potential_form_builder.Potential_Form_Builder` used to create potential instances.

  :returns: Potential callable that takes the product of a number of potential instances."""
  from atsim.potentials import product
  mod = _modifier_from_func_reduce("product", product, potential_forms, potential_form_builder)
  return mod

@modifier
def pow(potential_forms, potential_form_builder):
  """Modifier that takes two arguments and raises the first potential-form to the power of the second.

  e.g. pow(as.buck 1000.0 0.2 32.0, as.constant 2.0)

       Would square the given Buckingham potential.

  :param potential_forms: List of tuples that can be passed to `atsim.potentials.config._potential_form_builder.Potential_Form_Builder.create_potential_function()` to create potential callables.
  :param potential_form_builder: `atsim.potentials.config._potential_form_builder.Potential_Form_Builder` used to create potential instances.

  :returns: Potential callable that takes the product of a number of potential instances."""
  from atsim.potentials import pow
  mod = _modifier_from_func_reduce("pow", pow, potential_forms, potential_form_builder)
  return mod


class _Exp_Spline_Factory(object):
  """Helper class used by spline_potential() modifier to instantiate :class:`Exp_Spline` objects when `exp_spline` spline type is specified in config file."""


  spline_keyword = "exp_spline"

  def build_spline(self, detach_point, attach_point, spline_defn):
      if spline_defn.parameters:
        raise ConfigurationException("spline modifier 'exp_spline' middle potential form does not take any parameters. The following parameters were specified: {}".format(pot2.parameters))

      spline = Exp_Spline(detach_point, attach_point)
      return spline

class _Buck4_Spline_Factory(object):
  """Helper class used by spline_potential() modifier to instantiate :class:`Buck4_Spline` objects when `buck4_spline` spline type is specified in config file."""


  spline_keyword = "buck4_spline"

  def build_spline(self, detach_point, attach_point, spline_defn):
      if not len(spline_defn.parameters) == 1:
        raise ConfigurationException("spline modifier with 'buck4_spline' requires a single parameter to define r_min. The following parameters were specified: {}".format(pot2.parameters))

      r_min = spline_defn.parameters[0]

      if not r_min < attach_point.r and not r_min > detach_point.r:
        raise ConfigurationException("spline modifier with 'buck4_spline' r_min parameter does not lie between detach and attach values ({} < r_min < {}). r_min = {}".format(
          detach_point.r, attach_point.r, r_min))

      spline = Buck4_Spline(detach_point, attach_point, r_min)
      return spline



@modifier
def spline(potential_forms, potential_form_builder):
  """Modifier that smoothly splines between two potential forms by linking them with an intermediate spline.

  The `potential_forms` list must contain a single `PotentialFormInstanceTuple` entre. The tuple
  must only define three sections - potential form A -> spline -> potential form B. The extent of the spline is defined by where the middle region starts and potential form B starts.

  As a configuration string this might be defined as:

  `>0 as.zbl 14 8 >=0.8 exp_spline >=1.4 as.buck 180003 0.3 32.0`

  Which would create a `zbl` and Buckingham potential connected by a spline when `r` is between 0.8 and 1.4.

  At present the potential form label for the spline must be `exp_spline` or `buck4_spline`.

  :param potential_forms: List of tuples that can be passed to `atsim.potentials.config._potential_form_builder.Potential_Form_Builder.create_potential_function()` to create potential callables.
  :param potential_form_builder: `atsim.potentials.config._potential_form_builder.Potential_Form_Builder` used to create potential instances.

  :returns: Potential callable representing splined potentials."""

  spline_factories = [ _Exp_Spline_Factory(), _Buck4_Spline_Factory() ]

  if len(potential_forms) != 1:
    raise ConfigurationException("spline modifier only accepts a single multi range potential definition as its argument")

  logger = logging.getLogger(__name__).getChild("spline")

  # Extract the ranges and potential forms.
  pform = potential_forms[0]
  if pform.next is None:
    raise ConfigurationException("spline modifier requires three sub-potentials to be defined only one specified.")

  # ... set the 'next' attribute to None in pot1 and pot2
  pot1 = pform._replace(next = None)
  pot2 = pform.next._replace(next = None)

  allowed_spline_types = [s.spline_keyword for s in spline_factories]
  if not pot2.potential_form in allowed_spline_types:
    allowed_spline_types_str = ["'{}'".format(t) for t in allowed_spline_types]
    allowed_spline_types_str = ",".join(allowed_spline_types_str)
    raise ConfigurationException("spline modifier only accepts spline types {} for middle potential form. '{}' was found instead".format(
      allowed_spline_types_str,
      pot2.potential_form))

  if pform.next.next is None:
    raise ConfigurationException("spline modifier requires three sub-potentials to be defined only two specified.")

  if not pform.next.next.next is None:
    raise ConfigurationException("spline modifier requires three sub-potentials to be defined more than three have been given.")

  # ... we need to be able to evaluate the potential at the start point
  #     this allows values to be calculated at any point (as long as the potential's well behaved there)
  pot3_old_start = pform.next.next.start.start
  pot3 = pform.next.next._replace(
    start = MultiRangeDefinitionTuple(u">", float("-inf")),
    next = None)

  # Determine detachment point
  if not pot1.start.start < pot2.start.start:
    raise ConfigurationException("spline modifier range error. Start of 1st potential should be less than start of 2nd ! {} < {}".format(
      pot1.start.start,
      pot2.start.start))

  if not pot2.start.start < pot3_old_start:
    raise ConfigurationException("spline modifier range error. Start of 2nd potential (exp_spline) should be less than start of 3rd ! {} < {}".format(
      pot2.start.start,
      pot3_old_start))

  detach_point_r = pform.next.start.start
  attach_point_r = pform.next.next.start.start

  pot1_func = potential_form_builder.create_potential_function(pot1)
  pot3_func = potential_form_builder.create_potential_function(pot3)

  detach_point = Spline_Point(pot1_func, detach_point_r)
  attach_point = Spline_Point(pot3_func, attach_point_r)

  spline_factory = [s for s in spline_factories if s.spline_keyword == pot2.potential_form ][0]

  logger.debug("spline modifier: connecting '{}' with {} to '{}' in range {} to {}".format(
    pot1.potential_form,
    pot2.potential_form,
    pot2.potential_form,
    detach_point, attach_point))

  # Now build the spline object
  try:
    spline = spline_factory.build_spline(detach_point, attach_point, pot2)
  except ImportError as ie:
    raise_e = ConfigurationException("When using the '{}' spline type an additional package is required. Please install. {}".format(
        pot2.potential_form,
        str(ie.args[0])
    ))
    raise raise_e

  spot_obj = Custom_SplinePotential(spline)
  return spot_obj


@modifier
def trans(potential_forms, potential_form_builder):
  """Modifier that takes two arguments the first is a potential form and the second must be an
  instance of as.constant. This modifier performs the operation:

  potential_form(r+X)

  where X is the value parameter for the as.constant argument.

  e.g. trans(as.buck 1000.0 0.2 32.0, as.constant 2.0)

  :param potential_forms: List of tuples that can be passed to `atsim.potentials.config._potential_form_builder.Potential_Form_Builder.create_potential_function()` to create potential callables.
  :param potential_form_builder: `atsim.potentials.config._potential_form_builder.Potential_Form_Builder` used to create potential instances.

  :returns: Potential callable transforms the input to the potential-form."""
  if not len(potential_forms) == 2:
    raise ConfigurationException("trans() potential modifier only accepts two arguments")

  second_form = potential_forms[1]
  if second_form.potential_form != 'as.constant':
    raise ConfigurationException("the second argument to the trans() potential modifier must be 'as.constant' found {}".format(second_form))

  if len(second_form.parameters) != 1:
    raise ConfigurationException("the second parameter to trans(), 'as.constant' should have exactly one parameter defining shift. {} parameters found.".format(len(second_form.parameters)))

  logger = logging.getLogger(__name__).getChild("trans")


  potential_func = potential_form_builder.create_potential_function(potential_forms[0])
  trans_value = second_form.parameters[0]
  logger.debug("Creating trans() modifier with offset of {}".format(trans_value))
  def transformed(r):
    return potential_func(r+trans_value)

  if hasattr(potential_func, 'deriv'):
    def deriv(r):
      return potential_func.deriv(r+trans_value)
    transformed.deriv = deriv

  if hasattr(potential_func, 'deriv2'):
    def deriv2(r):
      return potential_func.deriv2(r+trans_value)
    transformed.deriv2 = deriv2

  return transformed
