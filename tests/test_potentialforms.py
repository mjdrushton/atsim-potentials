
from atsim.potentials import potentialfunctions, potentialforms, num_deriv

import pytest

import inspect

def test_polynomial():
  r = 1.2
  expect = 1.1 + 2.43*r + 1.82*r**2 - 2.4*r**3
  assert pytest.approx(potentialfunctions.polynomial(r, 1.1, 2.43, 1.82, -2.4)) == expect
  
  pf = potentialforms.polynomial(1.1, 2.43, 1.82, -2.4)
  assert pytest.approx(pf(r)) == expect

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

  # List of (args, expect) tuples
  expect = dict(
    bornmayer = ( (1000, 0.2), -1.6773131395 ),
    buck = ( (2000, 0.2, 32.0), 3.797931094),
    constant = ( (1.0,), 0),
    coul = ((-1, 2), -11.249609375),
    exp_spline = ((115.85848061468926, -554.8944421744643, 1103.3521708648423, -1086.8977963610287, 526.616672986912, -100.69329806728238, 0.0), -0.839792),
    hbnd = ((1000, 2000), 87.0415),
    lj = ((1.2,0.9), 0.534052),
    morse = ((0.1, 2.1, 1000), -10.7799643399247),
    polynomial = [ 
      ((1.0, 2.5), 2.5),
      ((1.0, 2.5, -3.0), -7.1),
      ((479.9553, -1372.5304, 1562.2233, -881.9685,  246.4347,  -27.2447), -2.10212479999882),
      (tuple(), 0.0),
    ],
    sqrt = ((3.0,), 1.1858541226),
    zbl = ((92, 8), -40.4739253056137),
    zero = (tuple(), 0.0)
  ) 

  for (potname, potargs) in expect.items():
    # potfunc = getattr(potentialfunctions, potname)

    if not isinstance(potargs, list):
      potargs = [potargs]

    for (args, e_val) in potargs:
      potfunc = caller.create_pot_func(potname, args)


      # Make sure we can still call the potential function
      caller.call_potfunc(potfunc, r, args)
      # potfunc(*r_args)

      actual = caller.call_potfunc(potfunc.deriv, r, args)
      assert pytest.approx(actual) == e_val

      # Check against the numerical derivative
      assert pytest.approx(actual, rel = 1e-2) == e_val


@pytest.mark.parametrize("module", [potentialfunctions])
def test_check_derivs(module):
  potfuncs = inspect.getmembers(module, potentialforms._iscallable)
  for (potname, potfunc) in potfuncs:
    if potname.startswith("_"):
      continue
    assert hasattr(potfunc, "deriv")