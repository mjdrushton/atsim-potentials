from __future__ import absolute_import
from __future__ import division
from builtins import object

import math

from ._util import gradient

def _splineExponential(startPoint, endPoint):
  """Returns a function of the form A + B*x + C*exp(x) + D*exp(-x) that
  passes through startPoint and endPoint with a gradient of startGradient
  at startPoint and endGradient at endPoint.

  :param startPoint: (x,y, dydx, ddydx) tuple where dydx is the gradient and ddydx is second derivative
  :param endPoint: (x,y, dydx, ddydx) tuple  where dydx is the gradient and ddydx is second derivative
  :return:  Tuple containing the interpolation function and the six polynomial coefficients"""

  #Calculate coefficients of the exponential function
  import numpy as np


  sx, sy, sdydx, sddydx = startPoint
  ex, ey, edydx, eddydx = endPoint

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

  coefficients = np.linalg.solve(A,B)

  B0 = float(coefficients[0])
  B1 = float(coefficients[1])
  B2 = float(coefficients[2])
  B3 = float( coefficients[3] )
  B4 = float( coefficients[4])
  B5 = float( coefficients[5])

  def potfunc(x):
    polynomial = B0 + B1*x + B2*x**2 + B3*x**3 + B4*x**4 + B5*x**5
    return math.exp(polynomial)

  return (potfunc, (B0, B1, B2, B3, B4, B5))


class SplinePotential(object):
  """Callable to allow splining of one potential to another"""

  def __init__(self, startPotential, endPotential, detachmentX, attachmentX):
    """Joins ``startPotential`` to ``endPotential`` using exponential spline of the form:

        .. math::

            U(r_{ij}) = \exp \left( B_0 + B_1 r_{ij} + B_2 r_{ij}^2 + B_3 r_{ij}^3 + B_4 r_{ij}^4 + B_5 r_{ij}^5 \\right)

    The spline coefficients :math:`B_{0...5}` can be obtained using the :meth:`.splineCoefficients` property.

    .. seealso::

      * :ref:`spline_interpolation`
      * :ref:`example_spline`


    :param startPotential: Function defining potential for rij <= detachmentX
    :param endPotential: Function defining potential for rij => attachmentX
    :param detachmentX: rij value at which startPotential should end
    :param attachmentX: rij value at which splines join endPotential"""


    startdydx = gradient(startPotential)
    startddydx = gradient(startdydx)

    enddydx = gradient(endPotential)
    endddydx = gradient(enddydx)

    self._startPotential = startPotential
    self._endPotential = endPotential
    self._attachmentX = attachmentX
    self._detachmentX = detachmentX

    detachmentY = startPotential(detachmentX)
    detachmentGradient = startdydx(detachmentX)
    detachmentGradientTwo = startddydx(detachmentX)


    attachmentY = endPotential(attachmentX)
    attachmentGradient = enddydx(attachmentX)
    attachmentGradient2 = endddydx(attachmentX)

    inter = 0.0

    #Can't take log of negative number, translate data upwards
    if detachmentY <=0.0 or attachmentY <= 0.0:
      inter = 1.0 - min([detachmentY, attachmentY])


    interfunc, self._splineCoefficients = _splineExponential(
      (detachmentX, detachmentY+inter, detachmentGradient, detachmentGradientTwo),
      (attachmentX, attachmentY+inter, attachmentGradient, attachmentGradient2))

    #Undo the effects of transposition
    def uninter(x):
      return interfunc(x) - inter

    self._interpolationFunction = uninter


  @property
  def startPotential(self):
    """:return:  Function defining potential for separations < ``detachmentX``"""
    return self._startPotential

  @property
  def endPotential(self):
    """:return:  Function defining potential for separations > ``attachmentX``"""
    return self._endPotential

  @property
  def interpolationFunction(self):
    """:return:  Exponential spline function connecting startPotential and endPotential for separations ``detachmentX`` < rij < ``attachmentX``"""
    return self._interpolationFunction

  @property
  def detachmentX(self):
    """:return:  Point at which spline should start"""
    return self._detachmentX

  @property
  def attachmentX(self):
    """:return:  Point at which spline should end"""
    return self._attachmentX

  @property
  def splineCoefficients(self):
    """:return:  Tuple containing the six coefficients of the spline polynomial"""
    return self._splineCoefficients

  def __call__(self, rij):
    """:param rij: separation at which to evaluate splined potential
       :return:  spline value"""
    if rij <= self.detachmentX:
      return self.startPotential(rij)
    elif rij >= self.attachmentX:
      return self.endPotential(rij)
    else:
      return self._interpolationFunction(rij)
