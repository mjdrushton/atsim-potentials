"""Tests for functions from the _util module"""

from atsim.potentials import gradient, deriv, num_deriv

import pytest


class GradientCallable(object):

  def __call__(self, r):
    return bare_function(r)

  def deriv(self, r):
    # If numerical derivative is being called thi should return 8.0 (r=1)
    # if this function is being called then 2.0
    return 2.0

  def deriv2(self, r):
    # If numerical derivative is taken then this should be  6.0
    # if this method is called then 1.0
    return 1.0

gradient_callable = GradientCallable()

def bare_function(r):
  return 3.0*r**2 + 2.0*r + 4.0
  

def test_num_deriv():
  expect = 8.0
  assert pytest.approx(num_deriv(1.0, gradient_callable)) == 8.0
  assert pytest.approx(num_deriv(1.0, bare_function)) == 8.0


def test_deriv():
  """Test the atsim.potentials._util.deriv function"""
  # Check that num_deriv is called when necessary
  assert pytest.approx(deriv(1.0, bare_function)) == 8.0

  # Check that deriv method is called where possible
  assert pytest.approx(deriv(1.0, gradient_callable )) == 2.0


def test_gradient():
  """Test the atsim.potentials._util.gradient function"""

  # Test with bare_function
  wrapped = gradient(bare_function)
  assert not hasattr(wrapped, "deriv")
  assert pytest.approx(wrapped(1.0)) == 8.0

  double_wrapped = gradient(wrapped)
  assert pytest.approx(double_wrapped(1.0), rel = 1e-4) == 6.0

  # Now test with objects that provide derivative methods.
  wrapped = gradient(gradient_callable)
  assert hasattr(wrapped, "deriv")
  assert wrapped(1.0) == 2.0

  double_wrapped = gradient(wrapped)
  assert not hasattr(double_wrapped, "deriv")
  assert double_wrapped(1.0) == 1.0