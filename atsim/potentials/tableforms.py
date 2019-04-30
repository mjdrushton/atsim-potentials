
class Cubic_Spline_Table_Form(object):
  """Potential form that takes tabulated data and returns interpolated values.

  This potential uses cubic spline interpolation. It is simply a wrapper around
  the scipy.interpolate.InterpolatedUnivariateSpline class"""
  
  """Identifier used in config files through `interpolation` directive"""
  config_label  = "cubic_spline"

  #Tag used during introspection of class so that it is automatically registered with Potential_Form_Registry
  is_potential = True

  def __init__(self, x_data, y_data):
    """Creates a potential form by linking x,y data points by cubic splines.

    :param x_data: List of x data values.
    :param y_data: List of y data values."""
    
    from scipy.interpolate import InterpolatedUnivariateSpline
    # ext =1 means that a value of zero is returned outside the data range. 
    self._interpolant = InterpolatedUnivariateSpline(x_data, y_data, ext =1)
    self._deriv = self._interpolant.derivative()
    self._deriv2 = self._deriv.derivative()

  @property
  def interpolant(self):
    """This class is a wrapper around instances of scipy.interpolate.InterpolatedUnivariateSpline 
    This property returns the scipy object used internally"""
    return self._interpolant

  def __call__(self, x):
    """:return: interpolated value at x"""
    return float(self._interpolant(x))


  def deriv(self, x):
    """:return: derivative of potential form at x"""
    return float(self._deriv(x))

  def deriv2(self, x):
    """:return: second derivative of potential form at x"""
    return float(self._deriv2(x))