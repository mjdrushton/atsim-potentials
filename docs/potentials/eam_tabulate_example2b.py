#! /usr/bin/env python

from atsim.potentials import writeTABEAM
from atsim.potentials import Potential
from atsim.potentials import EAMPotential

import math

def makeFunc(a, b, r_e, c):
  # Creates functions of the form used for density function.
  # Functional form also forms components of pair potential.
  def func(r):
    return (a * math.exp(-b*(r/r_e -1)))/(1+(r/r_e - c)**20.0)
  return func


def makePairPotAA(A, gamma, r_e, kappa,
  # Function factory that returns functions parameterised for homogeneous pair interactions
                  B, omega, lamda):
  f1 = makeFunc(A, gamma, r_e, kappa)
  f2 = makeFunc(B, omega, r_e, lamda)
  def func(r):
    return f1(r) - f2(r)
  return func


def makePairPotAB(dens_a, phi_aa, dens_b, phi_bb):
  # Function factory that returns functions parameterised for heterogeneous pair interactions
  def func(r):
    return 0.5 * ( (dens_b(r)/dens_a(r) * phi_aa(r)) + (dens_a(r)/dens_b(r) * phi_bb(r)) )
  return func

def makeEmbed(rho_e, rho_s, F_ni, F_i, F_e, eta):
  # Function factory returning parameterised embedding function.
  rho_n = 0.85*rho_e
  rho_0 = 1.15*rho_e

  def e1(rho):
    return sum([F_ni[i] * (rho/rho_n - 1)**float(i) for i in xrange(4)])

  def e2(rho):
    return sum([F_i[i] * (rho/rho_e - 1)**float(i) for i in xrange(4)])

  def e3(rho):
    return F_e * (1.0 - eta*math.log(rho/rho_s)) * (rho/rho_s)**eta

  def func(rho):
    if rho < rho_n:
      return e1(rho)
    elif rho_n <= rho < rho_0:
      return e2(rho)
    return e3(rho)
  return func

def makePotentialObjects():
  # Potential parameters
  r_eCu     = 2.556162
  f_eCu     = 1.554485
  gamma_Cu  = 8.127620
  omega_Cu  = 4.334731
  A_Cu      = 0.396620
  B_Cu      = 0.548085
  kappa_Cu  = 0.308782
  lambda_Cu = 0.756515

  rho_e_Cu  = 21.175871
  rho_s_Cu  = 21.175395
  F_ni_Cu   = [-2.170269, -0.263788, 1.088878, -0.817603]
  F_i_Cu    = [-2.19, 0.0, 0.561830, -2.100595]
  eta_Cu     = 0.310490
  F_e_Cu    = -2.186568

  r_eAl     = 2.863924
  f_eAl     = 1.403115
  gamma_Al  = 6.613165
  omega_Al  = 3.527021
  # A_Al      = 0.134873
  A_Al      = 0.314873
  B_Al      = 0.365551
  kappa_Al  = 0.379846
  lambda_Al = 0.759692

  rho_e_Al  =  20.418205
  rho_s_Al  =  23.195740
  F_ni_Al   = [-2.807602, -0.301435, 1.258562, -1.247604]
  F_i_Al    = [-2.83, 0.0, 0.622245, -2.488244]
  eta_Al     = 0.785902
  F_e_Al    = -2.824528

  # Define the density functions
  dens_Cu   = makeFunc(f_eCu, omega_Cu, r_eCu, lambda_Cu )
  dens_Al   = makeFunc(f_eAl, omega_Al, r_eAl,  lambda_Al )

  # Finally, define embedding functions for each species
  embed_Cu = makeEmbed(rho_e_Cu, rho_s_Cu, F_ni_Cu, F_i_Cu, F_e_Cu, eta_Cu)
  embed_Al = makeEmbed(rho_e_Al, rho_s_Al, F_ni_Al, F_i_Al, F_e_Al, eta_Al)

  # Wrap them in EAMPotential objects
  eamPotentials = [
    EAMPotential("Al", 13, 26.98, embed_Al, dens_Al),
    EAMPotential("Cu", 29, 63.55, embed_Cu, dens_Cu)]

  # Define pair functions
  pair_CuCu = makePairPotAA(A_Cu, gamma_Cu, r_eCu, kappa_Cu,
                           B_Cu, omega_Cu, lambda_Cu)

  pair_AlAl = makePairPotAA(A_Al, gamma_Al, r_eAl, kappa_Al,
                           B_Al, omega_Al, lambda_Al)

  pair_AlCu = makePairPotAB(dens_Cu, pair_CuCu, dens_Al, pair_AlAl)

  # Wrap them in Potential objects
  pairPotentials = [
      Potential('Al', 'Al', pair_AlAl),
      Potential('Cu', 'Cu', pair_CuCu),
      Potential('Al', 'Cu', pair_AlCu)]

  return eamPotentials, pairPotentials

def main():
  eamPotentials, pairPotentials = makePotentialObjects()

  # Perform tabulation
  # Make tabulation
  nrho = 2000
  drho = 0.05

  nr = 2000
  dr = 0.003

  with open("TABEAM", 'wb') as outfile:
    writeTABEAM(
      nrho, drho,
      nr, dr,
      eamPotentials,
      pairPotentials,
      out= outfile)

if __name__ == '__main__':
  main()
