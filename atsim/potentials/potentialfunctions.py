"""Functions for different potential forms."""

from __future__ import division

import math

def buck(r, A, rho, C):
  """Buckingham potential form

  .. math ::

    U(r_{ij}) = A \exp \left( \\frac{-r_{ij}}{\\rho} \\right) - \\frac{C}{r_{ij}^6}

  :param r: Atomic separation.
  :param A: Buckingham A parameter
  :param rho: Buckingham rho parameter :math:`\\rho`
  :param C: Buckingham C parameter
  :return: Potential energy."""
  return A * math.exp(-r/rho) - (C/r**6)

def bornmayer(r, A, rho):
  """Born-Mayer potential form

  .. math ::

    U(r_{ij}) = A \exp \left( \\frac{-r_{ij}}{\\rho} \\right)

  :param r: Atomic separation.
  :param A: Potential parameter
  :param rho: Potential parameter :math:`\\rho`
  :return: Potential energy"""
  return buck(r, A,rho,0.0)

def coul(r, qi, qj):
  """Coulomb potential (including :math:`4\pi \epsilon_0` term).

  .. math ::

    U(r_{ij}) = \\frac{ q_i q_j }{4 \pi \epsilon_0 r_{ij} }

  .. note:: Constant value appropriate for :math:`r_{ij}` in angstroms and energy in eV.

  :param r: Atomic separation.
  :param qi: Charge on species i
  :param qj: Charge on species j

  :return: Potential energy"""
  return (qi * qj)/(4.0*math.pi*0.0055264*r)

def constant(r, constant):
  """Function that simply returns the function that is passed to it.

  :param r: Separation.
  :param constant: The value to be returned by this function.

  :return: The value passed in as `constant`"""
  return constant

def hbnd(r, A,B):
  """DL_POLY hbnd form:

  .. math::

    U(r_{ij}) = \\frac{A}{r_{ij}^{12}} - \\frac{B}{r_{ij}^{10}}


  :param r: Atomic separation.
  :param A: Potential A parameter
  :param B: Potentials' B parameter
  :return: Potential energy"""
  return (A/r**12) - (B/r**10)

def lj(r, epsilon, sigma):
  """Lennard-Jones 12-6 potential.

  .. math ::

    U(r_{ij}) = 4 \epsilon \left( \\frac{\sigma^{12}}{r_{ij}^{12}} - \\frac{\sigma^6}{r_{ij}^6} \\right)

  :param r: Atomic separation.
  :param epsilon: Epsilon parameter :math:`\epsilon`
  :param sigma: Sigma parameter :math:`\sigma`
  :return: Potential energy"""
  return 4.0*epsilon*( (sigma**12/r**12) - (sigma**6/r**6))

def morse(r, gamma, r_star, D):
  """Morse potential parametrised with gamma, r_start and D.

  Potential function is (where r is interatomic separation):

  .. math ::

    U(r_{ij}) = D  \left[ \exp \left( -2 \gamma (r_{ij} - r_*) \\right) - 2 \\exp \left( -\gamma (r - r_*) \\right) \\right]


  :param r: Atomic separation.
  :param gamma: Potential parameter :math:`\gamma`
  :param r_star: Potential parameter :math:`r_*`
  :param D: Potential parameter
  :return: Potential energy"""
  return D*(math.exp(-2.0*gamma*(r-r_star)) - 2.0*math.exp(-gamma*(r-r_star)))


def sqrt(r,G):
  """Square root function.

  Potential function is:

  .. math ::

    U(r_{ij}) = G\sqrt{r_{ij}}

  :param r: Variable
  :param G: Parameter :math:`G`
  :return: G multiplied by the square root of r."""
  return G*math.sqrt(r)

def zbl(r, z1,z2):
  """ZBL potential.

  Ziegler-Biersack-Littmark screened nuclear repulsion for describing high energy interactions.

  :param r: Atomic separation.
  :param z1: Atomic number of species i
  :param z2: Atomic number of species j
  :return: Potential energy"""
  Ck1=0.1818
  Ck2=0.5099
  Ck3=0.2802
  Ck4=0.02817
  Bk1=3.2
  Bk2=0.9423
  Bk3=0.4029
  Bk4=0.2016
  a=(0.8854*0.529)/(float(z1)**0.23 + float(z2)**0.23)
  return 14.39942 * (z1*z2)/r * (Ck1*math.exp((-Bk1*r)/a) + Ck2*math.exp((-Bk2*r)/a) + Ck3*math.exp((-Bk3*r)/a) + Ck4*math.exp((-Bk4*r)/a))


def zero(r):
  """Function that returns 0.0 for any value of r.

  :param r: Separation.
  
  :return: Returns 0.0"""
  return 0.0

def exp_spline(r, B0, B1, B2, B3, B4, B5, C):
  """Exponential spline function (as used in splining routines).

    .. math::

            U(r_{ij}) = \exp \left( B_0 + B_1 r_{ij} + B_2 r_{ij}^2 + B_3 r_{ij}^3 + B_4 r_{ij}^4 + B_5 r_{ij}^5 \\right) + C

  :param B0: Spline coefficient
  :param B1: Spline coefficient
  :param B2: Spline coefficient
  :param B3: Spline coefficient
  :param B4: Spline coefficient
  :param B5: Spline coefficient
  :param C: C parameter
  :param r: Separation:

  :return: Splined values"""
  return math.exp(B0 + B1*r + B2*r**2 + B3*r**3 + B4*r**4 + B5*r**5) + C

