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
  return isinstance(obj, Callable)

class _FunctionFactory(object):

  def __init__(self, func):
    self._func = func

  def __call__(self, *args):
    return _rpartial(self._func, *args)

def _populate_module():
  potfuncs = inspect.getmembers(potentialfunctions, _iscallable)

  currmodule = sys.modules[__name__]
  for name, func in potfuncs:
    funcfactory = _FunctionFactory(func)
    setattr(currmodule, name, funcfactory)

_populate_module()
