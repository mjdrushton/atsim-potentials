# -*- coding: utf-8 -*-
import math
import os

from io import StringIO


import sys

from ._eam_potential import EAMPotential # noqa


def _writeHeader(outfile, nrho, drho, nr, dr, cutoff, title, atomicNumber, mass, latticeConstant, latticeType):
  print(title, file=outfile)
  print(u"%d %f %f %s" % (atomicNumber, mass, latticeConstant, latticeType), file=outfile)
  print(u"%d %f %d %f %f" % (nrho, drho, nr, dr, cutoff), file=outfile)

def _writeValueBlock(outfile, values):
  numbertemplate = u" % 20.16e"
  i = 0
  for value in values:
    i +=1
    if value == None :
      if i != 1:
        outfile.write(os.linesep)
      i = 0
      continue

    outfile.write(numbertemplate % value)

    if i % 5 == 0:
      i = 0
      outfile.write(os.linesep)
      continue

def _writeSetFLHeader(nrho, drho, nr, dr, cutoff, eampots, comments, out):
  #Line 1-3: comments
  workout = StringIO()

  newcomments = list(comments)
  newcomments.extend(['', '', ''])
  newcomments = newcomments[:3]

  print(os.linesep.join(newcomments), file=workout)

  #Line 4: ntypes
  ntypes = [u"%d" % len(eampots)]
  typestrings = [ eampot.species for eampot in eampots ]
  ntypes.extend(typestrings)
  ntypes = u" ".join(ntypes)

  print(ntypes, file=workout)

  #Line 5: nrho drho nr dr rcutoff
  numbertemplate = u" %20.16e"
  templ = u"%%d %s %%d %s %s" % ((numbertemplate,)*3)
  print(templ % (nrho, drho, nr, dr, cutoff), file=workout)

  #Dump everything into out
  out.write(workout.getvalue())

#Writes elementblock header
def _writeSetFLElementHeader(eampot, out):
  workout = StringIO()
  print(u"%d %20.16e %20.16e %s" % (eampot.atomicNumber, eampot.mass, eampot.latticeConstant, eampot.latticeType), file=workout)
  out.write(workout.getvalue())

#Writes element block embedding function
def _writeSetFLEmbeddingFunction(nrho, drho, eampot, out):
  workout = StringIO()
  for i in range(nrho):
    rho = float(i)*drho
    val = eampot.embeddingFunction(rho)
    print(u"% 20.16e" % val, file=workout)
  out.write(workout.getvalue())

def _writeDensityFunction(func, nr, dr, out):
  for i in range(nr):
    r = float(i)*dr
    val = func(r)
    print(u"% 20.16e" % val, file=out)

#Writes single density function to element block
def _writeSetFLDensityFunction(eampot, eampots, nr, dr, out):
  workout = StringIO()
  _writeDensityFunction(eampot.electronDensityFunction, nr, dr, workout)
  #Dump into out
  out.write(workout.getvalue())


#Writes multiple density functions (one per element) for each element block
def _writeSetFLDensityFunctionFinnisSinclair(eampot, eampots, nr, dr, out):
  workout = StringIO()

  for otherpot in eampots:
    densFunc = otherpot.electronDensityFunction[eampot.species]
    _writeDensityFunction(densFunc, nr, dr, workout)

  #Dump into out
  out.write(workout.getvalue())

def _writeSetFLPairPots(nr, dr, eampots, pairpots, out, scale_r = True):
  workout = StringIO()

  def pairkey(a,b):
    k = [a,b]
    k.sort()
    return tuple(k)

  class ZeroPair(object):
    def energy(self, rij):
      return 0.0
  zeroPair = ZeroPair()
  #Make a dictionary of available pair pots
  pairpotsdict = {}
  for pp in pairpots:
    pairpotsdict[pairkey(pp.speciesA,pp.speciesB)] = pp

  #Make list of required pots
  for i in range(len(eampots)):
    for j in range(i+1):
      k = pairkey(eampots[i].species, eampots[j].species)
      pp = pairpotsdict.get(k, zeroPair)

      for k in range(nr):
        r = float(k) * dr
        val = pp.energy(r)
        if scale_r:
          val *= r
        print(u"% 20.16e" % val, file=workout)
  out.write(workout.getvalue())

def _writeSetFL(
  nrho, drho,
  nr, dr,
  cutoff,
  eampots,
  pairpots,
  comments,
  out,
  writeDensityFunction):

  workout = StringIO()
  _writeSetFLHeader(nrho, drho, nr, dr, cutoff, eampots, comments, workout)

  #Write element block
  for eampot in eampots:
    eleblockout = StringIO()
    _writeSetFLElementHeader(eampot, eleblockout)
    _writeSetFLEmbeddingFunction(nrho, drho, eampot, eleblockout)
    writeDensityFunction(eampot, eampots, nr, dr, eleblockout)
    workout.write(eleblockout.getvalue())

  #Write pair potentials
  _writeSetFLPairPots(nr, dr, eampots, pairpots, workout)

  #Dump working output file into out
  out.write(workout.getvalue())

def writeFuncFL(
    nrho, drho,
    nr, dr,
    eampots,
    pairpots,
    out = sys.stdout,
    title = ""):
  """Creates a DYNAMO ``funcfl`` formatted file suitable for use with lammps `pair_style eam <http://lammps.sandia.gov/doc/pair_eam.html>`_
  potential form. For the `pair_style eam/alloy <http://lammps.sandia.gov/doc/pair_eam.html>`_ see :func:`~atsim.potentials.writeSetFL`.

  .. seealso ::

    For a working example using this function see :ref:`eam_example_1`


  :param nrho: Number of points used to describe embedding function
  :type nrho: int

  :param drho: Step size between rho values used to describe embedding function
  :type drho: float

  :param nr:   Number of points used for the pair-potential, and density functions
  :type nr: int

  :param dr:   Step size between r values in effective charge and density functions
  :type dr: float

  :param eampots: List containing a single :class:`~atsim.potentials.EAMPotential` instance for species to be tabulated.
  :type eampots: list

  :param pairpots: List containing a single :class:`~atsim.potentials.PairPotential` instance for the X-X interaction (where X is the species represented by EAMPotential in ``eampots`` list)
  :type pairpots: list

  :param out:  Python file object to which eam table file will be written
  :type out: file object

  :param title: Title to be written as table file header
  :type title: str"""

  cutoff = dr * (nr - 1)
  eampot = eampots[0]
  pairpot = pairpots[0]

  workingfile = StringIO()
  _writeHeader(workingfile,
               nrho, drho,
               nr, dr, cutoff,
               title,
               eampot.atomicNumber,
               eampot.mass,
               eampot.latticeConstant,
               eampot.latticeType)

  #Build a combined list of values
  rhos = [ float(x) * drho for x in range(nrho) ]
  embeds = [ eampot.embeddingFunction(rho) for rho in rhos]
  separations = [ float(x) * dr for x in range(nr) ]
  densities = [ eampot.electronDensityFunction(sep) for sep in separations]

  # Use pair potential to calculate energy * separation in ev Å
  charges = [ pairpot.energy(sep)*sep for sep in separations]
  # convert these into units of hartree * bohr radius
  charges = [ charge * 1.0/27.2 * 1.0/0.529 for charge in charges ]

  # TODO: If more precise conversion values used, answers diverge from that given by GULP hence why following line commented.
  # charges = [ charge * 1.0/27.2116 * 1.0/0.529177249 for charge in charges ]

  # to support mixing rules take the sqrt of these values
  charges = [ math.sqrt(charge) for charge in charges ]

  embeds.append(None)
  charges.append(None)

  valuelist = []
  valuelist.extend(embeds)
  valuelist.extend(charges)
  valuelist.extend(densities)
  _writeValueBlock(workingfile, valuelist)

  #Dump working file into out
  out.write(workingfile.getvalue())

def writeSetFL(
  nrho, drho,
  nr, dr,
  eampots,
  pairpots,
  out = sys.stdout,
  comments = ["", "", ""],
  cutoff = None):
  """Creates EAM potential in the DYNAMO ``setfl`` format. This format is suitable for
  use with the ``LAMMPS`` `pair_style eam/alloy <http://lammps.sandia.gov/doc/pair_eam.html>`_.

  .. seealso ::

    For a working example using this function see :ref:`eam_example_2a`


  :param nrho: Number of points used to describe embedding function
  :type nrho: int
  :param drho: Increment used when tabulating embedding function
  :type drho: float
  :param nr: Number of points used to describe density and pair potentials
  :type nr: int
  :param dr: Separation increment used when tabulating density function and pair potentials
  :type dr: float
  :param eampots: Instances of lammps.writeEAMTable.EAMPotential() which encapsulate information about each species
  :type eampots: list
  :param pairpots: Instance of potentials.Potential, these describe repulsive pair potential component of EAM potential
  :type pairpots: list
  :param out: Python file object into which EAM potential data should be written
  :type out: file object
  :param comments: List containing three strings, these form the header of the created file
  :type comments: list
  :param cutoff: Pair potential and density cutoff, if None then value of ``nr`` * ``dr`` is used.
  :type cutoff: float
  """

  if not cutoff:
    cutoff = nr*dr

  #Specialise _writeSetFL to use _writeSetFLDensityFunction to write a single density function
  _writeSetFL(
    nrho, drho,
    nr, dr,
    cutoff,
    eampots,
    pairpots,
    comments,
    out,
    _writeSetFLDensityFunction)

def writeSetFLFinnisSinclair(
  nrho, drho,
  nr, dr,
  eampots,
  pairpots,
  out = sys.stdout,
  comments = ["", "", ""],
  cutoff = None):
  """Creates Finnis-Sinclar EAM potential in the DYNAMO ``setfl`` format. The format should be used with the
  ``LAMMPS`` `eam/fs pair_style <http://lammps.sandia.gov/doc/pair_eam.html>`_.

  The :class:`~atsim.potentials.EAMPotential` instances within the ``eampots`` list are expected to provide individual density functions
  for each species pair in the species being tabulated. See :meth:`atsim.potentials.EAMPotential.__init__` for how these are specified
  to the :class:`atsim.potentials.EAMPotential` constructor.

  .. seealso ::

    For a working example using this function see :ref:`eam_example_3a`



  :param nrho: Number of points used to describe embedding function
  :type nrho: int
  :param drho: Increment used when tabulating embedding function
  :type drho: float
  :param nr: Number of points used to describe density and pair potentials
  :type nr: int
  :param dr: Separation increment used when tabulating density function and pair potentials
  :type dr: float
  :param eampots: Instances of lammps.writeEAMTable.EAMPotential() which encapsulate information about each species
  :type eampots: list
  :param pairpots: Instance of potentials.Potential, these describe repulsive pair potential component of EAM potential
  :type pairpots: list
  :param out: Python file object into which EAM potential data should be written
  :type out: file object
  :param comments: List containing three strings, these form the header of the created file
  :type comments: list
  :param cutoff: Pair potential and density cutoff. If None then value of ``nr`` * ``dr`` is used.
  :type cutoff: float"""

  if not cutoff:
    cutoff = nr * dr

  #Specialise _writeSetFL to use _writeSetFLDensityFunctionFinnisSinclar to write multiple density functions
  _writeSetFL(
    nrho, drho,
    nr, dr,
    cutoff,
    eampots,
    pairpots,
    comments,
    out,
    _writeSetFLDensityFunctionFinnisSinclair)
