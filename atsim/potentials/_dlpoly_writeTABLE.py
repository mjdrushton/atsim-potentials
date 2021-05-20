"""Module for creating DL_POLY TABLE files"""
import sys

from io import StringIO

class WritePotentialException(Exception):
  """Exception class that can be thrown by functions in this module"""
  pass


def _writePotential(potential, cutoff, gridPoints, meshResolution, out ):
  """Given a writeTABLE.Potential object, will write it to the given stream (out)
  in the correct DL_POLY TABLE file format.

  @param potential Potential object
  @param cutoff Short-range potential cutoff
  @param gridPoints Number of points used to describe potential
  @param meshResolution distance between grid points
  @param out Python object supporting write() to which output is sent"""

  #Check that number of grid points is divisible by 4
  if gridPoints%4 != 0:
    raise WritePotentialException("The number of rows in a DL_POLY TABLE file needs to be divisible by 4. Number of rows specified = {} ".format(gridPoints))

  outputbuilder = StringIO()

  #Generate output
  #Write the potential's header
  outputbuilder.write(u"%(atom1)8s%(atom2)8s\n" % { 'atom1' : potential.speciesA, 'atom2' : potential.speciesB})

  #Write the data records
  dataTemplate = u" % 14.7e"
  dataTemplate = dataTemplate+ u" % 14.7e" * 3
  dataTemplate = dataTemplate + u"\n"

  #First, do the energies
  l = []
  r=0.0
  for i in range(gridPoints):
    r += meshResolution
    l.append(potential.energy(r))

    if len(l) == 4:
      #List has 4 elements, dump a row
      outputbuilder.write(dataTemplate % tuple(l))
      #Reset the list
      l = []

  #Now, do the forces
  l = []
  r = 0.0
  for i in range(gridPoints):
    r += meshResolution
    l.append(_calculateForce(potential, r))

    if len(l) == 4:
      #List has 4 elements, dump a row
      outputbuilder.write(dataTemplate % tuple(l))
      #Reset the list
      l = []

  #Dump the output to out
  out.write(outputbuilder.getvalue())

def _calculateForce(pot, r):
  """Calls pot.force for separation (r) and returns DL_POLY -r dU/dr values rather than
  the dU/dr value normally returned by potentials.Potential.force() method

  @param pot potential from which force should be calculated
  @param r Separation at which force should be calculated

  @return -r dU/dr"""
  dUdr = pot.force(r)
  return r * dUdr

def _writeTableHeader(delpot, cutpot, ngrid, out):
  """Function responsible for creating header to DL_POLY TABLE file

  @param delpot Distance increment between grid points
  @param cutpot Potential short-range cutoff
  @param ngrid Number of grid points
  @param out Python stream object (supporting write()) to which output is sent"""
  outputbuilder = StringIO()
  outputbuilder.write(u" "*80 + u"\n")
  templParams = dict(delpot = delpot, cutpot = cutpot, ngrid = ngrid)
  outputbuilder.write(u"%(delpot)15.8e%(cutpot)15.8e%(ngrid)10d\n" % templParams)
  out.write(outputbuilder.getvalue())


def writePotentials(potentials, cutoff, gridPoints, out = sys.stdout):
  """Function used to convert a collection of Potential objects into a DL_POLY TABLE file.

  @param potentials Iterable containing Potential objects
  @param cutoff Short-range potential cutoff
  @param gridPoints Number of grid points used to tabulate potential
  @param out Python stream object (supporting write()) to which output is sent"""

  meshResolution = cutoff / (gridPoints-4.0)
  outputbuilder = StringIO()
  _writeTableHeader(meshResolution, cutoff, gridPoints, outputbuilder)

  for potential in potentials:
    _writePotential(potential, cutoff, gridPoints, meshResolution, outputbuilder)

  #Write the contents of outputbuilder to the final destination
  out.write(outputbuilder.getvalue())
