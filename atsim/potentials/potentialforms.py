"""Functions representing different potential forms.

The functions contained herein are function factories returning
a function that takes separation as its sole argument."""
from __future__ import division

import math

def buck(A, rho, C):
  """Returns a Buckingham potential function for a given set of A, rho and C parameters

  .. math ::

    U(r_{ij}) = A \exp \left( \\frac{-r_{ij}}{\\rho} \\right) - \\frac{C}{r_{ij}^6}

  :param A: Buckingham A parameter
  :param rho: Buckingham rho parameter :math:`\\rho`
  :param C: Buckingham C parameter
  :return: Function that will evaulate Buckingham potential for given A, rho and C"""

  def potential(r):
    return A * math.exp(-r/rho) - (C/r**6)
  return potential

def bornmayer(A, rho):
  """Return a Born-Mayer potential function for the given parameters

  .. math ::

    U(r_{ij}) = A \exp \left( \\frac{-r_{ij}}{\\rho} \\right)


  :param A: Potential parameter
  :param rho: Potential parameter :math:`\\rho`
  :return: Function that will evaluate Born-Mayer potential for given A and rho parameters"""
  return buck(A,rho,0.0)

def coul(qi, qj):
  """Coulomb potential (including :math:`4\pi \epsilon_0` term).

  .. math ::

    U(r_{ij}) = \\frac{ q_i q_j }{4 \pi \epsilon_0 r_{ij} }

  .. note:: Constant value appropriate for :math:`r_{ij}` in angstroms and energy in eV.

  :param qi: Charge on species i
  :param qj: Charge on species j

  :return: Coulomb callable"""
  def func(r):
    return (qi * qj)/(4.0*math.pi*0.0055264*r)
  return func


def hbnd(A,B):
  """Returns a function with the DL_POLY hbnd form:

  .. math::

    U(r_{ij}) = \\frac{A}{r_{ij}^{12}} - \\frac{B}{r_{ij}^{10}}


  :param A: Potential A parameter
  :param B: Potentials' B parameter
  :return: Function that will evaulate energy for potential of hbnd form for given A and B parameters"""

  def potential(r):
    return (A/r**12) - (B/r**10)
  return potential

def lj(epsilon, sigma):
  """Lennard-Jones 12-6 potential.

  .. math ::

    U(r_{ij}) = 4 \epsilon \left( \\frac{\sigma^{12}}{r_{ij}^{12}} - \\frac{\sigma^6}{r_{ij}^6} \\right)

  :param epsilon: Epsilon parameter :math:`\epsilon`
  :param sigma: Sigma parameter :math:`\sigma`
  :return: Lennard Jones callable"""
  def potential(r):
    return 4.0*epsilon*( (sigma**12/r**12) - (sigma**6/r**6))
  return potential

def morse(gamma, r_star, D):
  """Return morse function parametrised with gamma, r_start and D.

  Potential function is (where r is interatomic separation):

  .. math ::

    U(r_{ij}) = D  \left[ \exp \left( -2 \gamma (r_{ij} - r_*) \\right) - 2 \\exp \left( -\gamma (r - r_*) \\right) \\right]


  :param gamma: Potential parameter :math:`\gamma`
  :param r_star: Potential parameter :math:`r_*`
  :param D: Potential parameter
  :return: Morse potential callable"""
  def f(r):
    return D*(math.exp(-2.0*gamma*(r-r_star)) - 2.0*math.exp(-gamma*(r-r_star)))
  return f

def zbl(z1,z2):
  """ZBL potential.

  Ziegler-Biersack-Littmark screened nuclear repulsion for describing high energy interactions.

  :param z1: Atomic number of species i
  :param z2: Atomic number of species j
  :return: ZBL callable"""
  def zblpotential(r):
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
  return zblpotential

