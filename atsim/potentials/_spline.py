from __future__ import absolute_import
from __future__ import division
from builtins import object

import math

from ._util import gradient
from .potentialfunctions import exp_spline



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


class Exp_Spline(object):
  """Class for represention splines of the form: 

  .. math::

            U(r_{ij}) = \exp \left( B_0 + B_1 r_{ij} + B_2 r_{ij}^2 + B_3 r_{ij}^3 + B_4 r_{ij}^4 + B_5 r_{ij}^5 \\right) + C

    The spline coefficients :math:`B_{0...5}` and `C` can be obtained using the :meth:`.spline_coefficients` property.
  
  """
  
  def __init__(self, detach_point, attach_point):
    """Create a callable object that is a the exponential spline between `detach_point` and 
    `attach_point`.

    :param detach_point: Instance of Spline_Point giving start of spline.
    :param attach_point: Instance of Spline_Point giving end of spline."""

    self._detach_point = detach_point
    self._attach_point = attach_point

    self._init_spline_coefficients()

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
    self._coefficients = tuple(coefficients)

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
    B0,B1,B2,B3,B4,B5,C = self.spline_coefficients
    # polynomial = B0 + B1*r + B2*r**2 + B3*r**3 + B4*r**4 + B5*r**5
    # v = math.exp(polynomial) + C
    # return v
    return exp_spline(r, B0, B1, B2, B3, B4, B5, C)

class SplinePotential(object):
  """Callable to allow splining of one potential to another"""

  def __init__(self, startPotential, endPotential, detachmentX, attachmentX):
    """Joins ``startPotential`` to ``endPotential`` using spline of form:

      .. math::

            U(r_{ij}) = \exp \left( B_0 + B_1 r_{ij} + B_2 r_{ij}^2 + B_3 r_{ij}^3 + B_4 r_{ij}^4 + B_5 r_{ij}^5 \\right) + C

    The spline coefficients :math:`B_{0...5}` can be obtained using the :meth:`.SplineCoefficients` property.
    
    .. seealso::

      * :ref:`spline_interpolation`
      * :ref:`example_spline`
      * Custom_Spline_Potential for object to allow a different functional form for the spline.

    :param startPotential: Function defining potential for rij <= detachmentX
    :param endPotential: Function defining potential for rij => attachmentX
    :param detachmentX: rij value at which startPotential should end
    :param attachmentX: rij value at which splines join endPotential"""

    self._detach_point = Spline_Point(startPotential, detachmentX)
    self._attach_point = Spline_Point(endPotential, attachmentX)
    self._interpolationFunction = Exp_Spline(self._detach_point, self._attach_point)


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
    """:return:  Exponential spline function connecting startPotential and endPotential for separations ``detachmentX`` < rij < ``attachmentX``"""
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
