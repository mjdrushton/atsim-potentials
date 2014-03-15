"""Module containing functions for creating DL_POLY TABEAM files"""

import sys

from _common import Potential

try:
  import cStringIO as StringIO
except ImportError:
  import StringIO

def _tabulateFunction(outputfile, func, numpoints, step):
  outputbuilder = StringIO.StringIO()
  row = []
  for i in xrange(numpoints):
    row.append("%f" % (func( float(i) * step)))
    if len(row) == 4:
      print >>outputbuilder, " ".join(row)
      row = []
  if row:
    print >>outputbuilder, " ".join(row)
  outputfile.write(outputbuilder.getvalue())

def _writeEmbeddingFunction(eampotential, nrho, drho, outfile):
  outputbuilder = StringIO.StringIO()
  print >>outputbuilder, "embe %s %d 0.0 %f" % (eampotential.species, nrho, (nrho-1) * float(drho))

  _tabulateFunction(outputbuilder, eampotential.embeddingFunction, nrho, drho)
  outfile.write(outputbuilder.getvalue())


def _writeDensityFunction(speciesA, speciesB, electronDensityFunction, nr, dr, outfile):
  outputbuilder = StringIO.StringIO()
  if speciesA and speciesB:
    print >>outputbuilder, "dens %s %s %d 0.0 %f" % (speciesA, speciesB, nr, (nr-1) * float(dr))
  else:
    print >>outputbuilder, "dens %s %d 0.0 %f" % (speciesA, nr, (nr-1) * float(dr))
  _tabulateFunction(outputbuilder, electronDensityFunction, nr, dr)
  outfile.write(outputbuilder.getvalue())

def _writePairPotential(pairPotential, nr, dr, outfile):
  outputbuilder = StringIO.StringIO()
  print >>outputbuilder, "pair %s %s %d 0.0 %f" % (pairPotential.speciesA, pairPotential.speciesB, nr, ((nr-1) * float(dr)))

  def potentialCallable(r):
    return pairPotential.energy(r)

  _tabulateFunction(outputbuilder, potentialCallable, nr, dr)
  outfile.write(outputbuilder.getvalue())

def _writePairPotentials(eamPotentials, pairPotentials, nr, dr, outfile):
  pairs = set()
  for i in eamPotentials:
    for j in eamPotentials:
      k = tuple(sorted([i.species, j.species]))
      pairs.add(k)

  pairPotDict = {}
  for pot in pairPotentials:
    k = tuple(sorted([pot.speciesA, pot.speciesB]))
    pairPotDict[k] = pot

  def nullfunc(rij):
    return 0.0

  for k in sorted(pairs):
    if not pairPotDict.has_key(k):
      pot = Potential(k[0], k[1], nullfunc)
    else:
      pot = pairPotDict[k]
    _writePairPotential(pot, nr, dr, outfile)


def _writeTitle(title, out):
  title = "%s%s" % (title, 100*' ')[:100]
  print >>out, title


def _writeTABEAM_exceptDensity(nrho, drho, nr, dr, eamPotentials, pairPotentials, title, numpots, outputbuilder):
  _writeTitle(title, outputbuilder)


  print >>outputbuilder, "%d" % (numpots,)

  #Write the pair potentials
  _writePairPotentials(eamPotentials, pairPotentials, nr, dr, outputbuilder)

  #Write the embedding functions
  for eampotential in eamPotentials:
    _writeEmbeddingFunction(eampotential, nrho, drho, outputbuilder)


def writeTABEAM(nrho, drho, nr, dr, eampots, pairpots, out = sys.stdout, title = ""):
  """Create ``TABEAM`` file for use with the ``DL_POLY`` simulation code.

  .. seealso ::

    For a working example using this function see :ref:`eam_example_2b`


  :param nrho:  Number of entries in tabulated embedding functions
  :type nrho: int
  :param drho:  Step size between consecutive embedding function entries
  :type drho: float
  :param nr:  Number of entries in tabulated pair potentials and density functions
  :type nr: int
  :param dr:  Step size between entries in tabulated pair potentials and density functions
  :type dr: float
  :param eampots: Potentials List of potentials.EAMPotential objects
  :type eam: list
  :param pair: Potentials List of potentials.Potential objects
  :type pairpots: list
  :param out:  Python file object to which TABEAM data should be written
  :type out: file object
  :param title:  Title of TABEAM file
  :type title: str"""

  #Write the title
  outputbuilder = StringIO.StringIO()
  numpots = len(eampots)
  numpots = numpots*(numpots+5)/2
  _writeTABEAM_exceptDensity(nrho, drho, nr, dr, eampots, pairpots, title, numpots, outputbuilder)

  #Write the density functions
  for eampotential in eampots:
    _writeDensityFunction(eampotential.species, None, eampotential.electronDensityFunction, nr, dr, outputbuilder)

  out.write(outputbuilder.getvalue())


def writeTABEAMFinnisSinclair(nrho, drho, nr, dr, eampots, pairpots, out = sys.stdout, title = ""):
  """Create Exended EAM variant of DL_POLY ``TABEAM`` file.

  The :class:`.EAMPotential` instances within the ``eampots`` list are expected to provide individual density functions
  for each species pair in the species being tabulated. See :meth:`.EAMPotential.__init__` for how these are specified
  to the :class:`.EAMPotential` constructors.

  .. note :: The Extended EAM variant for which this function creates ``TABEAM`` files (i.e. metal potential type = eeam) is only supported in DL_POLY versions >= 4.05.


  .. seealso ::

    For a working example using this function see :ref:`eam_example_3b`



  :param nrho:  Number of entries in tabulated embedding functions
  :type nrho: int
  :param drho:  Step size between consecutive embedding function entries
  :type drho: float
  :param nr:  Number of entries in tabulated pair potentials and density functions
  :type nr: int
  :param dr:  Step size between entries in tabulated pair potentials and density functions
  :type dr: float
  :param eampots: Potentials List of :class:`.EAMPotential` objects
  :type eam: list
  :param pairpots: Potentials List of :class:`.Potential` objects
  :type pairpots: list
  :param out:  Python file object to which ``TABEAM`` data should be written
  :type out: file object
  :param title:  Title of TABEAM file
  :type title: str"""

  #Write the title
  outputbuilder = StringIO.StringIO()
  numpots = len(eampots)
  numpots = 3*numpots*(numpots+1)/2
  _writeTABEAM_exceptDensity(nrho, drho, nr, dr, eampots, pairpots, title, numpots, outputbuilder)

  # Write the density functions
  speciesList = sorted([ep.species for ep in eampots])
  for eamPotential in eampots:
    speciesA = eamPotential.species
    for speciesB in speciesList:
      try:
        densityFunction = eamPotential.electronDensityFunction[speciesB]
      except KeyError:
        raise KeyError("Density function for '%s-%s' pair not specified in '%s' EAMPotential, density dictionary." % (speciesA, speciesB, speciesA))
      _writeDensityFunction(speciesA, speciesB, densityFunction, nr, dr, outputbuilder)
  out.write(outputbuilder.getvalue())
