"""A collection of classes and functions related to defining potentials"""
from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

from builtins import range

from builtins import object

from . import _tablereaders

from ._common import * # noqa
from ._spline import SplinePotential # noqa
from ._util import gradient # noqa
from ._dlpoly_writeTABEAM import writeTABEAM # noqa
from ._dlpoly_writeTABEAM import writeTABEAMFinnisSinclair # noqa

from ._lammpsWriteEAM import writeFuncFL # noqa
from ._lammpsWriteEAM import writeSetFL # noqa
from ._lammpsWriteEAM import writeSetFLFinnisSinclair # noqa

from .potentialforms import * # noqa

import sys


def plus(a,b):
  """Takes two functions and returns a third which when evaluated returns the result of a(r) + b(r)

  This function is useful for combining existing potentials.

  **Example:**

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

  Data is written to a text file as two columns (r and E) separated by spaces
  with no header.

  :param fileobj: Python file object into which data should be plotted
  :param lowx: X-axis lower value
  :param highx: X-axis upper value
  :param func: Function to be plotted
  :param steps: Number of data points to be plotted"""
  step = (highx - lowx) / float(steps)

  for i in range(steps):
    v = lowx + float(i)*step
    y = func(v)

    fileobj.write("{0} {1}\n".format(v,y))


def plot(filename, lowx, highx, func, steps=10000):
  """Convenience function for plotting the potential functions contained herein.

  Data is written to a text file as two columns (r and E) separated by spaces
  with no header.

  :param filename: File into which data should be plotted
  :param lowx: X-axis lower value
  :param highx: X-axis upper value
  :param func: Function to be plotted
  :param steps: Number of data points to be plotted"""

  with open(filename, 'w') as outfile:
    plotToFile(outfile, lowx, highx, func, steps)


def plotPotentialObject(filename, lowx, highx, potentialObject, steps=10000):
  """Convenience function for plotting energy of pair interactions
  given by instances of :class:`atsim.potentials.Potential` obtained by calling
  `potential` `.energy()` method.

  Data is written to a text file as two columns (r and E) separated by spaces
  with no header.

  :param filename: File into which data should be plotted
  :param lowx: X-axis lower value
  :param highx: X-axis upper value
  :param func: :class:`atsim.potentials.Potential` object.
  :param steps: Number of data points to be plotted"""

  with open(filename, 'w') as outfile:
    plotPotentialObjectToFile(outfile, lowx, highx, potentialObject, steps)


def plotPotentialObjectToFile(fileobj, lowx, highx, potentialObject, steps=10000):
  """Convenience function for plotting energy of pair interactions
  given by instances of :class:`atsim.potentials.Potential` obtained by calling
  `potential` `.energy()` method.

  Data is written to a text file as two columns (r and E) separated by spaces
  with no header.

  :param fileobj: Python file object into which data should be plotted
  :param lowx: X-axis lower value
  :param highx: X-axis upper value
  :param func: :class:`atsim.potentials.Potential` object.
  :param steps: Number of data points to be plotted"""

  def f(r):
    return potentialObject.energy(r)
  plotToFile(fileobj, lowx, highx, f, steps)


def _LAMMPS_writePotentials(potentialList, cutoff, gridPoints, out):
  """Wrapper function that adapts lammps.writeTABLE.writePotentials() to the API
  expected by potentials.writePotentials"""
  minr = cutoff / float(gridPoints)
  from ._lammps_writeTABLE import writePotentials
  writePotentials(potentialList,minr, cutoff, gridPoints, out)

class UnsupportedTabulationType(Exception):
  """Exception thrown by writePotentials() when unknown tabulation type specified"""
  pass

def writePotentials(outputType, potentialList, cutoff, gridPoints, out = sys.stdout):
  """Tabulates pair-potentials in formats suitable for multiple simulation codes.

  * The ``outputType`` parameter can be one of the following:

    * ``DL_POLY``:

      * This function creates output that can be written to a ``TABLE`` and used
        within DL_POLY.
      * for a working example see :ref:`Quick-Start: DL_POLY <quick_start>`.

    * ``LAMMPS``:

      * Creates files readable by LAMMPS `pair_style table <http://lammps.sandia.gov/doc/pair_table.html>`_
      * Each file can contain multiple potentials:

        * the block representing each potential has a title formed from the ``speciesA``
          and ``speciesB`` attributes of the :class:`~atsim.potentials.Potential`
          instance represented by the block. These are sorted into their natural
          order and separated by a hyphen to form the title.
        *  **Example:**

            * For a :class:`~atsim.potentials.Potential` where ``speciesA`` = Xe
              and ``speciesB`` = O the block title would be: ``O-Xe``.
            * If ``speciesA`` = B
              and ``speciesB`` = O the block title would be: ``B-O``.

        * within LAMMPSthe block title is used as the ``keyword`` argument to the
          `pair_style table <http://lammps.sandia.gov/doc/pair_table.html>`_
          ``pair_coeff`` directive.
      * For a working example see :ref:`quick_start_lammps`


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

  from ._dlpoly_writeTABLE import writePotentials as DLPOLY_writePotentials
  supportedTabulations = {
    'DL_POLY' : DLPOLY_writePotentials,
    'LAMMPS'  : _LAMMPS_writePotentials
    }

  if not outputType in supportedTabulations:
    raise UnsupportedTabulationType("Unsupported tabulation type: '%s' should be one of: %s" % (outputType, ",".join(sorted(supportedTabulations.keys()))))

  supportedTabulations[outputType](potentialList, cutoff, gridPoints, out)



