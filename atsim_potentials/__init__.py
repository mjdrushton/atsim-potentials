"""A collection of classes and functions related to defining potentials"""
import contextlib

import _tablereaders

from _common import * # noqa
from _spline import SplinePotential # noqa
from _util import gradient # noqa
from _dlpoly_writeTABEAM import writeTABEAM # noqa
from _dlpoly_writeTABEAM import writeTABEAMFinnisSinclair # noqa

from _lammpsWriteEAM import writeFuncFL # noqa
from _lammpsWriteEAM import writeSetFL # noqa
from _lammpsWriteEAM import writeSetFLFinnisSinclair # noqa

from potentialforms import * # noqa

import sys


def plus(a,b):
  """Takes two functions and returns a third which when evaluated returns the result of a(r) + b(r)

  This function is useful for combining existing potentials.

  Example::


    To combine :func:`.buck` and :func:`.hbnd` functions from the :mod:`.potentialsforms` module to give:

      .. code:: python

        A*(-r/rho) + C/r**6 + D/r**12 - E/r**10

    this function can then be used as follows:

      .. code:: python

        plus(buck(A,rho,C), hbnd(D,E))

  :param a: First callable
  :param b: Second callable
  :return: Function that when evaulated returns ``a(r) + b(r)``"""

  def potential(r):
    return a(r) + b(r)
  return potential


class TableReader(object):
  """Callable that allows pretabulated data to be used with a Potential object."""

  def __init__(self, fileobject):
    """The file passed to tablereader() is assumed to have one separation, energy pair per line, separated by
    a space.

    :param fileobject: Python file object from which separation, energy pairs should be read"""
    self._tablereader = _tablereaders.DatReader(fileobject)

  @property
  def datReader(self):
    """:return: _tablereaders.DatReader associated with this callable"""
    return self._tablereader

  def __call__(self, separation):
    return self._tablereader.getValue(separation)


def plotToFile(fileobj,lowx, highx, func, steps=10000):
  """Convenience function for plotting the potential functions contained herein.

  Data is written to a text file as two columns (r and E) separated by spaces.

  :param fileobj: Python file object into which data should be plotted
  :param lowx: X-axis lower value
  :param highx: X-axis upper value
  :param func: Function to be plotted
  :param steps: Number of data points to be plotted"""
  steps = float(steps)
  step = (highx - lowx)/steps

  for i in xrange(int(steps)):
    v = lowx + float(i)*step
    y = func(v)

    print >>fileobj, "%f %f" % (v, y)


def plot(filename, lowx, highx, func, steps=10000):
  """Convenience function for plotting the potential functions contained herein.

  Data is written to a text file as two columns (r and E) separated by spaces.

  :param filename: File into which data should be plotted
  :param lowx: X-axis lower value
  :param highx: X-axis upper value
  :param func: Function to be plotted
  :param steps: Number of data points to be plotted"""

  with contextlib.closing(open(filename, 'wb')) as outfile:
    plotToFile(outfile, lowx, highx, func, steps)


def _LAMMPS_writePotentials(potentialList, cutoff, gridPoints, out):
  """Wrapper function that adapts lammps.writeTABLE.writePotentials() to the API
  expected by potentials.writePotentials"""
  minr = cutoff/float(gridPoints)
  from _lammps_writeTABLE import writePotentials
  writePotentials(potentialList,minr, cutoff, gridPoints, out)

class UnsupportedTabulationType(Exception):
  """Exception thrown by writePotentials() when unknown tabulation type specified"""
  pass

def writePotentials(outputType, potentialList, cutoff, gridPoints, out = sys.stdout):
  """Tabulates pair-potentials in formats suitable for multiple simulation codes.

  :param outputType: The type of output that should be created can be one of: ``DL_POLY`` or ``LAMMPS``
  :type outputType: str

  :param potentialList: List of Potential objects to be tabulated.
  :type potentialList: list

  :param cutoff: Largest separation to be tabulated.
  :type cutoff: float

  :param gridPoints: Number of rows in tabulation.
  :type gridPoints: int

  :param out: Python file like object to which tabulation should be written
  :type out: file"""

  from _dlpoly_writeTABLE import writePotentials as DLPOLY_writePotentials
  supportedTabulations = {
    'DL_POLY' : DLPOLY_writePotentials,
    'LAMMPS'  : _LAMMPS_writePotentials
    }

  if not outputType in supportedTabulations:
    raise UnsupportedTabulationType("Unsupported tabulation type: '%s' should be one of: %s" % (outputType, ",".join(sorted(supportedTabulations.keys()))))

  supportedTabulations[outputType](potentialList, cutoff, gridPoints, out)



