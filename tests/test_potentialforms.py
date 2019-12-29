
from atsim.potentials import potentialfunctions, potentialforms, num_deriv, plus, gradient

import pytest

import inspect
import math
import types

def testBuck():
  r = 0.5
  a = 1388.773
  rho = 2.76
  c = 175.0

  potfunc = potentialforms.buck(a,rho,c)
  pytest.approx(-10041.34343) == potfunc(r)

def testExponential():
  r = 0.5
  A = 14.7
  n= 4

  potfunc = potentialforms.exponential(A, n)
  pytest.approx(0.91875) == potfunc(r)

def testHbnd():
  a = 300.0
  b = 20.0

  potfunc = potentialforms.hbnd(a,b)
  pytest.approx(1208320) == potfunc(0.5)

def testPlus():
  potfunc = plus( potentialforms.buck(1388.773, 2.76, 175.0), potentialforms.hbnd(300.0, 20.0))
  pytest.approx(1198278.65656831) == potfunc(0.5)

def testPlusDeriv():
  """Test that the plus() function returns wrapped functions that make use of analytical derivatives where available"""
  
  class AnalyticalDeriv(object):

    def __init__(self, m, deriv):
      self._m = m
      self._deriv = deriv

    def __call__(self, r):
      return self._m * r + 3.0

    def deriv(self, r):
      return self._deriv

  # First test that if neither function in the plus has an analytical derivative,
  # that the wrapped function also doesn't have a .deriv() method
  def f(r):
    return -4.0*r

  def g(r):
    return 10.0*r

  wrapped = plus(f,g)
  assert not hasattr(wrapped, "deriv")
  assert pytest.approx(wrapped(2.0)) == 12.0
  assert pytest.approx(6.0) == gradient(wrapped)(2.0)

  # Now test that if both functions provide analytical derivatives they are 
  # used and the wrapped function provides a deriv function.
  f_a = AnalyticalDeriv(-4.0, -24.0)
  g_a = AnalyticalDeriv(10.0, 100.0)

  wrapped = plus(f_a, g_a)
  assert hasattr(wrapped, "deriv")
  assert pytest.approx(wrapped(2.0)) == 18.0
  assert pytest.approx(6.0) == num_deriv(2.0, wrapped)
  assert pytest.approx(76.0) == wrapped.deriv(2.0)
  assert pytest.approx(76.0) == gradient(wrapped)(2.0)

  # Now test that if one of the functions doesn't provide a deriv - the
  # one that does is still used. 
  wrapped = plus(f_a, g)
  assert hasattr(wrapped, "deriv")
  assert pytest.approx(wrapped(2.0)) == 15.0
  assert pytest.approx(6.0) == num_deriv(2.0, wrapped)
  assert pytest.approx(-14.0) == wrapped.deriv(2.0)
  assert pytest.approx(-14.0) == gradient(wrapped)(2.0)

  # Now do tests for deriv2() method
  wrapped = plus(f,g)
  assert not hasattr(wrapped, "deriv2")

  r = 1.6

  A = 1000.0
  rho = 0.1
  C = 32

  def b_f(r):
    return A * math.exp(-r/rho) - C/r**6

  A_h = 100.0
  B_h = 50.0

  def h_f(r):
     return  (A_h/r**12) - (B_h/r**10)

  # Create functions with analytical derivatives 
  b_a = potentialforms.buck(A, rho, C)
  h_a = potentialforms.hbnd(A_h, B_h)

  assert hasattr(b_a , 'deriv2')
  assert hasattr(h_a,  'deriv2')
  assert not hasattr(b_f , 'deriv2')
  assert not hasattr(h_f,  'deriv2')

  expect_energy = pytest.approx(h_f(r) + b_f(r))

  wrapped_f = plus(h_f, b_f)
  assert not hasattr(wrapped_f, 'deriv2')
  assert expect_energy == wrapped_f(r)

  wrapped_a = plus(h_a, b_a)
  assert hasattr(wrapped_a, 'deriv2')
  assert expect_energy == wrapped_a(r)

  wrapped_m = plus(h_f, b_a)
  assert hasattr(wrapped_m, 'deriv2')
  assert expect_energy == wrapped_m(r)

  expect = pytest.approx(-29.1717612428203)
  expect_loose = pytest.approx(-29.1717612428203, rel = 1e-4)
  assert expect_loose == gradient(gradient(wrapped_f))(r)
  assert expect == gradient(gradient(wrapped_a))(r)
  assert expect_loose == gradient(gradient(wrapped_m))(r)

  assert expect == wrapped_a.deriv2(r)
  assert expect_loose == wrapped_m.deriv2(r)
  
  # Now let's deliberately change the deriv2 returned by h_a and b_a to make sure the analytical forms are being used.
  def d2(self, r):
    return 20.0

  def d2_2(self, r):
    return 30.0

  h_a.deriv2 = types.MethodType(d2, h_a)#, h_a.__class__)
  b_a.deriv2 = types.MethodType(d2_2, b_a)#, b_a.__class__)

  wrapped = plus(h_a, b_a)
  assert pytest.approx(20.0+30.0) == wrapped.deriv2(r)
  assert pytest.approx(20.0+30.0) == gradient(gradient(wrapped))(r)

def test_buck4():

  with pytest.raises(ImportError):
    from atsim.potentials.potentialfunctions import buck4

  from atsim.potentials.potentialforms import buck4

  A = 11272.6
  rho = 0.1363
  C = 134.0

  r_detach = 1.2
  r_min = 2.1
  r_attach = 2.6

  pot = buck4(A, rho, C, r_detach, r_min, r_attach)

  assert pytest.approx(0.8227238) == pot(1.3)
  assert pytest.approx(-0.662829) == pot(2.4)

  assert pytest.approx(potentialfunctions.bornmayer.deriv(1.0, A, rho)) == pot.deriv(1.0)
  assert pytest.approx(potentialfunctions.buck.deriv(3.0, 0.0, 0.1, C)) == pot.deriv(3.0)

  assert pytest.approx(potentialfunctions.bornmayer.deriv2(1.0, A, rho)) == pot.deriv2(1.0)
  assert pytest.approx(potentialfunctions.buck.deriv2(3.0, 0.0, 0.1, C)) == pot.deriv2(3.0)


def test_polynomial():
  r = 1.2
  expect = 1.1 + 2.43*r + 1.82*r**2 - 2.4*r**3
  assert pytest.approx(potentialfunctions.polynomial(r, 1.1, 2.43, 1.82, -2.4)) == expect
  
  pf = potentialforms.polynomial(1.1, 2.43, 1.82, -2.4)
  assert pytest.approx(pf(r)) == expect


def test_tang_toennies():
  # Old implementation
  def f2n(x, n):
    v = 0.0
    for k in range(2*n+1):
      v += x**float(k) / math.factorial(k)
    return 1.0 - math.exp(-x) * v


  def _ttSecondTerm(r, b, C_6, C_8, C_10):
    n_3 = f2n(b*r, 3) * C_6/r**6.0
    n_4 = f2n(b*r, 4) * C_8/r**8.0
    n_5 = f2n(b*r, 5) * C_10/r**10.0
    return n_3 + n_4 + n_5


  def tangToennies(A, b, C_6, C_8, C_10):
    def f(r):
      # Convert r from angstrom to a.u.
      f = A * math.exp(-b*r)
      s = _ttSecondTerm(r, b, C_6, C_8, C_10)
      v = f-s
      # Convert energy in a.u. to eV
      return v
    return f

  def auConvertWrap(func):
    def wrapped(r):
      r /= 0.5292
      v = func(r)
      return v * 27.211
    return wrapped

  A = 832.4
  b = 1.865
  C_6 = 129.6
  C_8 = 4187.0
  C_10 = 155500.0
  Kr_Kr_expect = auConvertWrap(tangToennies(A, b, C_6, C_8, C_10))
  d1 = gradient(Kr_Kr_expect)
  d2 = gradient(d1)

  # Check the sympy implementation
  sp_expression = potentialfunctions.tang_toennies._as_sympy()
  for r in range(2, 100):
    r = float(r) /10.0
    expect = pytest.approx(Kr_Kr_expect(r))
    actual = sp_expression.subs({"A" : A, "b": b, 
      "C_6" : C_6, "C_8" : C_8,
      "C_10" : C_10, "r" : r })
    assert expect == actual

  # New implementation
  for r in range(2, 100):
    r = float(r) /10.0

    expect = pytest.approx(Kr_Kr_expect(r))
    actual = potentialfunctions.tang_toennies(r, A, b, C_6, C_8, C_10)

    assert expect == actual

  # First and second derivatives
  for r in range(10, 40):
    r = float(r) /10.0

    expect = pytest.approx(d1(r), rel=1e-3)
    actual = potentialfunctions.tang_toennies.deriv(r, A, b, C_6, C_8, C_10)
    assert expect == actual

    expect = pytest.approx(d2(r), rel = 1e-3)
    actual = potentialfunctions.tang_toennies.deriv2(r, A, b, C_6, C_8, C_10)
    assert expect == actual


class PotentialFunctionsCaller(object):
  """Helper used in test_deriv, this is for testing the functions from atsim.potentialfunctions"""

  def create_pot_func(self, potname, args):
    potfunc = getattr(potentialfunctions, potname)
    return potfunc

  def call_potfunc(self, potfunc, r, args):
    r_args = [r]
    if args:
      r_args.extend(args)
    return potfunc(*r_args)

class PotentialFormsCaller(object):
  """Helper used in test_deriv, this is for testing the functions from atsim.potentialforms"""

  def create_pot_func(self, potname, args):
    potcls = getattr(potentialforms, potname)
    potfunc = potcls(*args)
    return potfunc

  def call_potfunc(self, potfunc, r, args):
    return potfunc(r)

@pytest.mark.parametrize("caller", [
  PotentialFunctionsCaller(),
  PotentialFormsCaller()])
def test_deriv(caller):
  """Check that potential functions have deriv() attributes"""

  r = 1.6

  # List of (args, deriv_expect, deriv2_expect) tuples
  expect = dict(
    bornmayer = ( (1000, 0.2), -1.6773131395, 8.38656569756279),
    buck = ( (2000, 0.2, 32.0), 3.797931094, -14.5193071119545),
    constant = ( (1.0,), 0, 0),
    coul = ((-1, 2), 11.249609375, -14.062024927003026),
    exp_spline = ((115.85848061468926, -554.8944421744643, 1103.3521708648423, -1086.8977963610287, 526.616672986912, -100.69329806728238, 0.0), -0.839792, 30.2521673632487),
    hbnd = ((1000, 2000), 87.0415, -565.103519534204),
    lj = ((1.2,0.9), 0.534052, -2.20102077321182),
    morse = ((0.1, 2.1, 1000), -10.7799643399247, 23.1814147955054),
    polynomial = [ 
      ((1.0, 2.5), 2.5, 0.0),
      ((1.0, 2.5, -3.0), -7.1, -6.0),
      ((479.9553, -1372.5304, 1562.2233, -881.9685,  246.4347,  -27.2447), -2.10212479999882, -3.86283999999978),
      (tuple(), 0.0, 0.0),
    ],
    tang_toennies = ((832.4, 1.865, 129.6, 4187.0, 155500.0), -271.564487468955, 973.169664930065),
    sqrt = ((3.0,), 1.1858541226, -0.370579413300982),
    zbl = ((92, 8), -40.4739253056137, 142.687518894378),
    zero = (tuple(), 0.0, 0.0),
    exponential = [
      ((1.0, 2), 3.2, 2),
      ((-5.8, 4), -95.0272, -178.176),
      ((3, 0.5), 1.1858541226, -0.370579413300982),
      ((1.2, -1), -0.46875, 0.5859375),
      ((1.2, -2), -0.5859375, 1.0986328125)
    ]
  ) 

  for (potname, potargs) in expect.items():
    # potfunc = getattr(potentialfunctions, potname)

    if not isinstance(potargs, list):
      potargs = [potargs]

    for (args, e_deriv, e_deriv2) in potargs:
      potfunc = caller.create_pot_func(potname, args)

      # Make sure we can still call the potential function
      caller.call_potfunc(potfunc, r, args)
      # potfunc(*r_args)

      # First derivative
      actual = caller.call_potfunc(potfunc.deriv, r, args)
      assert pytest.approx(actual) == e_deriv

      # Check against the numerical derivative
      assert pytest.approx(num_deriv(r, lambda r: caller.call_potfunc(potfunc, r, args)), rel = 1e-2) == e_deriv

      # Second derivative
      actual = caller.call_potfunc(potfunc.deriv2, r, args)
      assert pytest.approx(actual) == e_deriv2

      # Check against the numerical derivative
      assert pytest.approx(num_deriv(r, lambda r: caller.call_potfunc(potfunc.deriv, r, args)), rel = 1e-2) == e_deriv2


@pytest.mark.parametrize("module", [potentialfunctions])
def test_check_derivs(module):
  potfuncs = inspect.getmembers(module, potentialforms._iscallable)
  for (potname, potfunc) in potfuncs:
    if potname.startswith("_"):
      continue
    assert hasattr(potfunc, "deriv")
    assert hasattr(potfunc, "deriv2")