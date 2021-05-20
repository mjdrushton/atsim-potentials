from io import StringIO

import os
import sys

def _writeSinglePotential(pot, minr, maxr, gridPoints, out):
  """Create a lammps tabulated potential for the given potentials.Potential
  object.

  The potential name within the LAMMPs file is of the form pot.speciesA-potspeciesB

  @param pot potentials.Potential instance
  @param minr Starting separation that should be used for potential
  @param maxr Maximum separation for which potential should be written
  @param gridPoints number of points to be included in tabulated form of potential
  @param out Python file object supporting write() method to which tabulated
             potential will be written"""
  #Write the section header
  sbuild = StringIO()

  print(u"%s-%s" % (pot.speciesA, pot.speciesB), file=sbuild)
  print(u"N %(gridpoints)d R %(minr).8f %(maxr).8f" % { 'gridpoints' : gridPoints,
                                                                 'minr' : minr,
                                                                 'maxr' : maxr }, file=sbuild)
  print(u"", file=sbuild)
  #Write the body of the potential
  for n in range(1,gridPoints+1):
    r = minr + float(n-1)* (maxr - minr) / (float(gridPoints) -1)
    energy = pot.energy(r)
    force = pot.force(r)

    print(u"%(n)s %(r).8f %(energy).8f %(force).8f" % { 'n' : n,
                                                                  'r' : r,
                                                                  'energy' :  energy,
                                                                  'force' : force }, file=sbuild)
  out.write(sbuild.getvalue())

def writePotentials(potentials, minr, maxr, gridPoints, out = sys.stdout):
  """Formats potentials.Potentials instances into a LAMMPs tabulated potential.

  @param potentials Iterable containing potentials.Potential instances
  @param minr Separation at which tabulation of potential should start
  @param maxr Maximum separation for tabulation of potential
  @param gridPoints number of grid points between minr and maxr at which each potential should be sampled
  @param out Python file object into which tabulated potentials should be written"""

  potlines = []
  for potential in potentials:
    sbuild = StringIO()
    _writeSinglePotential(potential, minr, maxr, gridPoints, sbuild)
    potlines.append(sbuild.getvalue())
  out.write(os.linesep.join(potlines))

