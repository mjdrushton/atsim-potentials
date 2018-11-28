"""Functions representing different potential forms.

The functions contained herein are function factories returning
a function that takes separation as its sole argument."""

from . import potentialfunctions
import inspect
import sys

from ._util import _rpartial

try:
  from collections.abc import Callable
except ImportError:
  from collections import Callable


def _iscallable(obj):
  return (not obj is Callable) and isinstance(obj, Callable) and (not inspect.isclass(obj))

class _FunctionFactory(object):

  def __init__(self, func):
    self._func = func

  def __call__(self, *args):
    wrapper = _rpartial(self._func, *args)
    if hasattr(self._func, "deriv"):
      wrapper.deriv = _rpartial(self._func.deriv, *args)
    return wrapper


def _populate_module():
  potfuncs = inspect.getmembers(potentialfunctions, _iscallable)

  currmodule = sys.modules[__name__]
  for name, func in potfuncs:
    funcfactory = _FunctionFactory(func)
    setattr(currmodule, name, funcfactory)

_populate_module()
