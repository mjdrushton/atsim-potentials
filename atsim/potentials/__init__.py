"""A collection of classes and functions related to defining potentials"""

from . import _tablereaders

from ._potential import Potential
from ._eam_potential import EAMPotential

from .spline import SplinePotential # noqa
from ._util import gradient, num_deriv, deriv # noqa
from ._dlpoly_writeTABEAM import writeTABEAM # noqa
from ._dlpoly_writeTABEAM import writeTABEAMFinnisSinclair # noqa

from ._lammpsWriteEAM import writeFuncFL # noqa
from ._lammpsWriteEAM import writeSetFL # noqa
from ._lammpsWriteEAM import writeSetFLFinnisSinclair # noqa

from .potentialforms import * # noqa
from . import potentialfunctions
from . import spline
from . import pair_tabulation
from . import eam_tabulation
from . import tableforms

from ._multi_range_potential_form import create_Multi_Range_Potential_Form, Multi_Range_Defn

import sys
import math

def plus(a,b):
  """Takes two functions and returns a third which when evaluated returns the result of a(r) + b(r)

  This function is useful for combining existing potentials.

  **Derivatives:**

  If either of the potential callables (`a` and `b`) provide a .deriv() method the function returned by
  `plus()` will also have a `.deriv()` method. This allows analytical derivatives to be specified. If
  only one of `a` or `b` provide `.deriv()` then the derivative of the other callable will be evaluated
  numerically.

  If neither function has a .deriv() method then the function returned here will also *not* have a .deriv()
  method.

  **Example:**

    To combine :func:`~atsim.potentials.potentialforms.buck` and :func:`~atsim.potentials.potentialforms.hbnd` functions from the :mod:`atsim.potentials.potentialforms` module to give:

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

  # Set derivatives
  if hasattr(a, 'deriv') or hasattr(b, 'deriv'):
    deriv_a = gradient(a)
    deriv_b = gradient(b)
    def deriv(r):
      return deriv_a(r) + deriv_b(r)
    potential.deriv = deriv

    if hasattr(deriv_a, 'deriv') or hasattr(deriv_b, 'deriv'):
      deriv2_a = gradient(deriv_a)
      deriv2_b = gradient(deriv_b)
      def deriv2(r):
        return deriv2_a(r) + deriv2_b(r)
      potential.deriv2 = deriv2

  return potential

def product(a,b):
  """Takes two callables and returns a third which when evaluated returns the result of a(r) * b(r)

  This function is useful for combining existing potentials.

  **Derivatives:**

  If either of the potential callables (`a` and `b`) provide a .deriv() method the function returned by
  `product()` will also have a `.deriv()` method. This allows analytical derivatives to be specified. If
  only one of `a` or `b` provide `.deriv()` then the derivative of the other callable will be evaluated
  numerically.

  If neither function has a .deriv() method then the function returned here will also *not* have a .deriv()
  method.

  :param a: First callable
  :param b: Second callable
  :return: Function that when evaulated returns ``a(r) * b(r)``"""

  def potential(r):
    return a(r) * b(r)

  # Set derivatives
  if hasattr(a, 'deriv') or hasattr(b, 'deriv'):
    deriv_a = gradient(a)
    deriv_b = gradient(b)
    def deriv(r):
      return a(r) * deriv_b(r) + b(r) * deriv_a(r)
    potential.deriv = deriv

    if hasattr(deriv_a, 'deriv') or hasattr(deriv_b, 'deriv'):
      deriv2_a = gradient(deriv_a)
      deriv2_b = gradient(deriv_b)
      def deriv2(r):
        return deriv2_a(r) * b(r) + 2.0*deriv_a(r)*deriv_b(r) + a(r)*deriv2_b(r)
      potential.deriv2 = deriv2
  return potential

def pow(a,b):
  """Takes two callables and returns a third which when evaluated returns the result of a(r)**b(r)

  This function is useful for combining existing potentials.

  **Derivatives:**

  If either of the potential callables (`a` and `b`) provide a .deriv() method the function returned by
  `pow()` will also have a `.deriv()` method. This allows analytical derivatives to be specified. If
  only one of `a` or `b` provide `.deriv()` then the derivative of the other callable will be evaluated
  numerically.

  If neither function has a .deriv() method then the function returned here will also *not* have a .deriv()
  method.

  :param a: First callable
  :param b: Second callable
  :return: Function that when evaulated returns ``a(r)**b(r)`` (a to the power of b)"""

  def potential(r):
    return a(r)**b(r)

  # Set derivatives
  if hasattr(a, 'deriv') or hasattr(b, 'deriv'):
    deriv_a = gradient(a)
    deriv_b = gradient(b)
    def deriv(r):
      ar = a(r)
      return potential(r) * (deriv_b(r) * math.log(ar) + b(r) * deriv_a(r)/ar)
    potential.deriv = deriv

    if hasattr(deriv_a, 'deriv') or hasattr(deriv_b, 'deriv'):
      deriv2_a = gradient(deriv_a)
      deriv2_b = gradient(deriv_b)
      def deriv2(r):
        ar = a(r)
        br = b(r)
        dr = deriv(r)
        p = potential(r)
        da = deriv_a(r)
        db = deriv_b(r)
        d2a = deriv2_a(r)
        d2b = deriv2_b(r)

        # value = (deriv_b(r)*log(a(r)) + b(r)*deriv_a(r)/a(r))*deriv(r) + (math.log(a(r))*deriv2_b(r) + b(r)*deriv2_a(r)/a(r) + deriv_a(r)*deriv2_b(r)/a(r) + deriv_b(r)*deriv2_a(r)/a(r) - b(r)*deriv_a(r)*deriv2_a(r)/a(r)**2)*potential(r)
        value = (db*math.log(ar) + (br*da)/ar)*dr + (math.log(ar)*d2b + (br*d2a)/ar + (da*db)/ar + (db*da)/ar - (br*da*da)/(ar**2))*p
        return value
      potential.deriv2 = deriv2
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
      * for a working example see :ref:`Quick-Start: DL_POLY <quick-start>`.

    * ``GULP``:

      * Creates output for the `GULP code <https://nanochemistry.curtin.edu.au/gulp/>`_
      * Output is in the form of a series of `spline potential forms <https://nanochemistry.curtin.edu.au/gulp/help/new_help_40_txt.html#spline>`_
      * The generated file can be loaded into GULP using the `library command <https://nanochemistry.curtin.edu.au/gulp/help/new_help_40_txt.html#library>`_

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

  from .pair_tabulation import DLPoly_PairTabulation, LAMMPS_PairTabulation, GULP_PairTabulation

  supportedTabulations = {
    'DL_POLY' : DLPoly_PairTabulation,
    'LAMMPS'  : LAMMPS_PairTabulation,
    'GULP'    : GULP_PairTabulation
    }

  if not outputType in supportedTabulations:
    raise UnsupportedTabulationType("Unsupported tabulation type: '{}' should be one of: {}".format(
      outputType, ",".join(sorted(supportedTabulations.keys()))
    ))

  tabulation = supportedTabulations[outputType](potentialList,cutoff, gridPoints)
  tabulation.write(out)
  



