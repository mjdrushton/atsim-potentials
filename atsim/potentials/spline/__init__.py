import math

from .._util import gradient

from ..potentialforms import polynomial, exp_spline


class Spline_Point(object):
  """Class for the attachment and detachment points of potential objects and region to be splined"""

  def __init__(self, potential_function, r):
    """Create object which evaluates the value and 1st and 2nd derivatives of `potential` at a point `r`.

    :param potential_function: Function to be evaluated, this is a callable accepting a single argument.
    :param r: The point at which function should be evaluated"""
    self._potential_function = potential_function
    self._r = r

    self._deriv_callable = gradient(self._potential_function)
    self._deriv2_callable = gradient(self._deriv_callable)

  @property
  def potential_function(self):
    """Potential function"""
    return self._potential_function

  @property
  def r(self):
    """Value at which splining takes place"""
    return self._r

  @property
  def v(self):
    """Value of `potential_function` at `r`"""
    return self.potential_function(self.r)

  @property
  def deriv(self):
    """First derivative of `potential_function`: dv/dr(r)"""
    return self._deriv_callable(self.r)

  @property
  def deriv2(self):
    """Second derivative of `potential_function`: d2v/dr^2(r)"""
    return self._deriv2_callable(self.r)
  
  @property
  def deriv_callable(self):
    return self._deriv_callable

  @property
  def deriv2_callable(self):
    return self._deriv2_callable


class Exp_Spline(object):
  """Class for represention splines of the form: 

  .. math::

            U(r_{ij}) = \exp \left( B_0 + B_1 r_{ij} + B_2 r_{ij}^2 + B_3 r_{ij}^3 + B_4 r_{ij}^4 + B_5 r_{ij}^5 \\right) + C

  
  The spline coefficients :math:`B_{0...5}` and `C` can be obtained using the :meth:`~atsim.potentials.Exp_Spline.spline_coefficients` property.
  
  """
  
  def __init__(self, detach_point, attach_point):
    """Create a callable object that is an exponential spline between `detach_point` and 
    `attach_point`.

    :param detach_point: Instance of Spline_Point giving start of spline.
    :param attach_point: Instance of Spline_Point giving end of spline."""

    self._detach_point = detach_point
    self._attach_point = attach_point

    self._coefficients = self._init_spline_coefficients()
    self._spline_callable = exp_spline(*self.spline_coefficients)

  def _init_spline_coefficients(self):
    #Calculate coefficients of the exponential function
    import numpy as np

    sx = self.detach_point.r
    sy = self.detach_point.v
    sdydx = self.detach_point.deriv
    sddydx = self.detach_point.deriv2

    ex = self.attach_point.r
    ey = self.attach_point.v
    edydx = self.attach_point.deriv
    eddydx = self.attach_point.deriv2

    inter = 0.0

    #Can't take log of negative number, translate data upwards
    if sy <=0.0 or ey <= 0.0:
      inter = 1.0 - min([sy, ey])
      sy += inter
      ey += inter
      inter = -inter

    A = np.array(
                [[1.0 , sx  , sx**2  , sx**3     , sx**4      , sx**5 ]       ,
                [1.0  , ex  , ex**2  , ex**3     , ex**4      , ex**5]        ,
                [0.0  , 1.0 , 2.0*sx , 3.0*sx**2 , 4.0*sx**3  , 5.0*sx**4]    ,
                [0.0  , 1.0 , 2.0*ex , 3.0*ex**2 , 4.0*ex**3  , 5.0*ex**4]    ,
                [0.0  , 0.0 , 2.0    , 6.0*sx    , 12.0*sx**2 , 20.0*sx**3]   ,
                [0.0  , 0.0 , 2.0    , 6.0*ex    , 12.0*ex**2 , 20.0*ex**3]])

    B = np.array([
        math.log(sy),
        math.log(ey),
        sdydx / sy,
        edydx / ey,
        (sddydx/sy)-((sdydx**2)/(sy**2)),
        (eddydx/ey)-((edydx**2)/(ey**2))
      ])

    coefficients = [float(c) for c in np.linalg.solve(A,B)]
    coefficients.append(inter)
    return tuple(coefficients)

  @property
  def detach_point(self):
    """Spline_Point giving start of splined region"""
    return self._detach_point

  @property
  def attach_point(self):
    """Spline_Point giving end of splined region"""
    return self._attach_point

  @property
  def spline_coefficients(self):
    """Coefficients for spline_function"""
    return self._coefficients

  def __call__(self, r):
    return self._spline_callable(r)

  def deriv(self, r):
    return self._spline_callable.deriv(r)

  def deriv2(self, r):
    return self._spline_callable.deriv2(r)

class Buck4_Spline(object):
  """Class for representing the splined part of the four ranged Buckingham potential.

  Between the detachment point and `r_min` this is a 5th order polynomial:

  .. math::

    U(r_{ij}) = A_0 + A_1 r_{ij} + A_2 r_{ij}^2 + A_3 r_{ij}^3 + A_4 r_{ij}^4 + A_5 r_{ij}^5

  and between `r_min` and the re-attachment point a 3rd order spline is used:

  .. math::

    U(r_{ij}) = B0 + B_1 r_{ij} + B_2 r_{ij}^2 + B_3 r_{ij}^3

  The spline coefficients :math:`A_{0..5}` and :math:`B_{0..3}` are solved such that the 
  the spline values match with the potential functions at the detach and re-attachment points and r_min.
  They are continuous in their first and second derivatives across these points and where the two
  splines meet at `r_min`. Finally, the derivative at `r_min` is set to be 0 with the aim of creating a
  minimum."""


  def __init__(self, detach_point, attach_point, r_min):
    """Create a callable object that represents the Buckingham-4 type spline between `detach_point` and 
    `attach_point`.

    :param detach_point: Instance of :class:`Spline_Point` giving start of spline.
    :param attach_point: Instance of :class:`Spline_Point` giving end of spline.
    :param r_min: Minimum value to be formed between the detach and reattachment points."""

    self._detach_point = detach_point
    self._attach_point = attach_point
    self._r_min = r_min

    self._init_spline_coefficients()


  def _init_spline_coefficients(self):
    r_dp = self.detach_point.r
    r_dp2 = r_dp**2
    r_dp3 = r_dp**3
    r_dp4 = r_dp**4
    r_dp5 = r_dp**5

    r_min = self.r_min
    r_min2 = r_min**2
    r_min3 = r_min**3
    r_min4 = r_min**4
    r_min5 = r_min**5

    r_ap = self.attach_point.r
    r_ap2 = r_ap**2
    r_ap3 = r_ap**3

    import numpy as np

    M = [ 
      1, r_dp, r_dp2  , r_dp3   , r_dp4    , r_dp5    , 0 , 0     , 0       , 0        ,
      0, 1   , 2*r_dp , 3*r_dp2 , 4*r_dp3  , 5*r_dp4  , 0 , 0     , 0       , 0        ,
      0, 0   , 2      , 6*r_dp  , 12*r_dp2 , 20*r_dp3 , 0 , 0     , 0       , 0        ,
      
      0, 1   , 2*r_min, 3*r_min2, 4*r_min3 , 5*r_min4 , 0 ,0      ,0        ,0         ,
      1,r_min, r_min2 , r_min3  , r_min4   , r_min5   , -1, -r_min, -r_min2 , -r_min3  ,
      0,1    , 2*r_min, 3*r_min2, 4*r_min3 , 5*r_min4 , 0 , -1    , -2*r_min, -3*r_min2,
      0, 0   , 2      , 6*r_min , 12*r_min2, 20*r_min3, 0 , 0     , -2      , -6*r_min ,

      0, 0   , 0      , 0       , 0        , 0        , 1 , r_ap  , r_ap2   , r_ap3    ,
      0, 0   , 0      , 0       , 0        , 0        , 0 , 1     , 2*r_ap  , 3*r_ap2  ,
      0, 0   , 0      , 0       , 0        , 0        , 0 , 0     , 2       , 6*r_ap]

    M = np.reshape(M, (10,10))

    V = [
      self.detach_point.v,
      self.detach_point.deriv,
      self.detach_point.deriv2,
      0,
      0,
      0,
      0,
      self.attach_point.v,
      self.attach_point.deriv,
      self.attach_point.deriv2 ]

    V = np.reshape(V, (10,1))

    coefficients = np.linalg.solve(M,V)
    coefficients = coefficients.flatten().tolist()

    self._spline5 = polynomial(*coefficients[:6])
    self._spline3 = polynomial(*coefficients[6:])


  @property
  def detach_point(self):
    """Spline_Point giving start of splined region"""
    return self._detach_point

  @property
  def attach_point(self):
    """Spline_Point giving end of splined region"""
    return self._attach_point

  @property
  def r_min(self):
    """Position of minimum"""
    return self._r_min

  @property
  def spline_coefficients(self):
    """Spline coefficients as list of form [A_0, A_1, A_2, A_3, A_4, A_5, B_0, B_1, B_2, B_3]"""
    return self.spline5.args + self.spline3.args

  @property
  def spline5(self):
    """Callable (atsim.potentials.potentialfunctions.polynomial) object representing the fifth order section of the buck4 spline - between `detach_point` and `r_min`"""
    return self._spline5

  @property
  def spline3(self):
    """Callable (atsim.potentials.potentialfunctions.polynomial) object representing the fifth order section of the buck4 spline - between `detach_point` and `r_min`"""
    return self._spline3

  def _which_spline(self, r):
    if r < self.r_min:
      return self.spline5
    else:
      return self.spline3

  def __call__(self, r):
    spline = self._which_spline(r)
    return spline(r)

  def deriv(self, r):
    return self._which_spline(r).deriv(r)

  def deriv2(self, r):
    return self._which_spline(r).deriv2(r)


class Custom_SplinePotential(object):
  """Callable to allow splining of one potential to another"""

  def __init__(self, spline):
    """Adapts spline objects such as :class:`Exp_Spline` and :class:`Buck4_Spline` to have the same interface as :class:`SplinePotential`
      
    :param spline: Instance of a spline object such as :class:`Exp_Spline` and :class:`Buck4_Spline`

    """

    self._spline = self._interpolationFunction = spline
    self._detach_point = self._spline.detach_point
    self._attach_point = self._spline.attach_point

    self._init_deriv()

  def _init_deriv(self):
    # This isn't technically used as a point - but Spline_Point is used to conveniently set up the gradient callables.
    self._inter_point = Spline_Point(self._interpolationFunction, None)
    # Expose a .deriv() method if any of components provide one natively.
    if hasattr(self._detach_point.potential_function, "deriv") or hasattr(self._attach_point.potential_function, "deriv") or hasattr(self._inter_point.potential_function, "deriv"):
      def deriv(self, r):
        return self._deriv(r)
      self.deriv = deriv.__get__(self)

    # Expose a .deriv2() method if any of components provide one natively.
    if hasattr(self._detach_point.potential_function, "deriv2") or hasattr(self._attach_point.potential_function, "deriv2") or hasattr(self._inter_point.potential_function, "deriv2"):
      def deriv2(self, r):
        return self._deriv2(r)
      self.deriv2 = deriv2.__get__(self)


  # The _deriv and _deriv2 methods are made public as deriv() and deriv2() by _init_deriv() if any of the component functions provide a deriv or deriv2 method.
  def _deriv(self, rij):
    if rij <= self.detachmentX:
      return self._detach_point.deriv_callable(rij)
    elif rij >= self.attachmentX:
      return self._attach_point.deriv_callable(rij)
    else:
      return self._inter_point.deriv_callable(rij)

  def _deriv2(self, rij):
    if rij <= self.detachmentX:
      return self._detach_point.deriv2_callable(rij)
    elif rij >= self.attachmentX:
      return self._attach_point.deriv2_callable(rij)
    else:
      return self._inter_point.deriv2_callable(rij)

  @property
  def startPotential(self):
    """:return:  Function defining potential for separations < ``detachmentX``"""
    return self._detach_point.potential_function

  @property
  def endPotential(self):
    """:return:  Function defining potential for separations > ``attachmentX``"""
    return self._attach_point.potential_function

  @property
  def interpolationFunction(self):
    """:return:  Spline object connecting startPotential and endPotential for separations ``detachmentX`` < rij < ``attachmentX``"""
    return self._interpolationFunction

  @property
  def detachmentX(self):
    """:return:  Point at which spline should start"""
    return self._detach_point.r

  @property
  def attachmentX(self):
    """:return:  Point at which spline should end"""
    return self._attach_point.r

  @property
  def splineCoefficients(self):
    """:return:  Tuple containing the seven coefficients of the spline polynomial"""
    return self._interpolationFunction.spline_coefficients

  def __call__(self, rij):
    """:param rij: separation at which to evaluate splined potential
       :return:  spline value"""
    if rij <= self.detachmentX:
      return self.startPotential(rij)
    elif rij >= self.attachmentX:
      return self.endPotential(rij)
    else:
      return self._interpolationFunction(rij)

class SplinePotential(Custom_SplinePotential):
  """Callable to allow splining of one potential to another using an exponential spline"""

  def __init__(self, startPotential, endPotential, detachmentX, attachmentX):
    """Joins ``startPotential`` to ``endPotential`` using spline of form:

      .. math::

            U(r_{ij}) = \exp \left( B_0 + B_1 r_{ij} + B_2 r_{ij}^2 + B_3 r_{ij}^3 + B_4 r_{ij}^4 + B_5 r_{ij}^5 \\right) + C

    The spline coefficients :math:`B_{0...5}` can be obtained using the :meth:`~atsim.potnetials.SplinePotential.splineCoefficients` property.
    
    .. seealso::

      * :ref:`spline_interpolation`
      * :ref:`example_spline`
      * :class:`Exp_Spline` - for details of the class which actually performs splining.
      
    :param startPotential: Function defining potential for rij <= detachmentX
    :param endPotential: Function defining potential for rij => attachmentX
    :param detachmentX: rij value at which startPotential should end
    :param attachmentX: rij value at which splines join endPotential"""

    detach_point = Spline_Point(startPotential, detachmentX)
    attach_point = Spline_Point(endPotential, attachmentX)
    spline = Exp_Spline(detach_point, attach_point)

    super(SplinePotential, self).__init__(spline)


class Buck4_SplinePotential(Custom_SplinePotential):
  """Callable to allow splining of one potential to another using the Buck4 spline type"""

  def __init__(self, startPotential, endPotential, detachmentX, attachmentX, r_min):
    """Joins ``startPotential`` to ``endPotential`` using :class:`Buck4_Spline`. This class is provided for
    convenience to make instantiating a :class:`Custom_SplinePotential` using this spline type more straightforward.

    .. seealso::

      * :ref:`spline_interpolation`
      * :class:`Buck4_Spline` - for details of the class which actually performs splining.
      
    :param startPotential: Function defining potential for rij <= detachmentX
    :param endPotential: Function defining potential for rij => attachmentX
    :param detachmentX: rij value at which startPotential should end
    :param attachmentX: rij value at which splines join endPotential
    :param r_min: rij value for Buck4 spline minimum."""

    detach_point = Spline_Point(startPotential, detachmentX)
    attach_point = Spline_Point(endPotential, attachmentX)
    spline = Buck4_Spline(detach_point, attach_point, r_min)

    super(Buck4_SplinePotential, self).__init__(spline)
