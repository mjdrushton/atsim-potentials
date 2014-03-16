import unittest
import StringIO
import contextlib
import os
import math
import subprocess
import shutil


from atsim.potentials import Potential, EAMPotential, buck, TableReader
from atsim import potentials

from _tempfiletestcase import TempfileTestCase
from _runlammps import needsLAMMPS, extractLAMMPSEnergy, runLAMMPS

def _getResourceDirectory():
  """Returns path to resources used by this test module (currently assumed to be sub-directory
  of test module called resources)"""
  return os.path.join(os.path.dirname(__file__), 'lammps_resources')

def _openResource(resourcename):
  """Returns file object for file in resources subdirectory"""
  return open( os.path.join(_getResourceDirectory(), resourcename), 'r')

def _floatiterator(fobj):
  for line in fobj:
    line = line[:-1]
    tokens = line.split()
    for token in tokens:
      yield float(token)

def _parseEAMTable(infile):
  """Parses a DYNAMO funcfl table file into a dictionary.

  @param infile Python file object containing eam file information
  @return A dictionary containing parsed tokens from the table"""

  outDict = {}
  outDict['title'] = infile.next()[:-1]

  line = infile.next()[:-1]
  atomicNumber, mass, latticeConstant, latticeType = line.split()
  outDict['atomicNumber'] = int(atomicNumber)
  outDict['mass'] = float(mass)
  outDict['latticeConstant'] = float(latticeConstant)
  outDict['latticeType'] = latticeType

  line = infile.next()[:-1]
  nrho, drho, nr, dr, cutoff = line.split()
  outDict['nRho'] = nrho = int(nrho)
  outDict['dRho'] = drho = float(drho)
  outDict['nr'] = nr = int(nr)
  outDict['dr'] = dr = float(dr)
  outDict['cutoff'] = cutoff = float(cutoff)

  #Read the embedding function values
  iter = _floatiterator(infile)
  for i in xrange(nrho):
    outDict.setdefault('embeddingFunction', []).append(iter.next())

  #Read the effective charge function values
  for i in xrange(nr):
    outDict.setdefault('effectiveChargeFunction', []).append(iter.next())

  #Read the density function values
  for i in xrange(nr):
    outDict.setdefault('densityFunction', []).append(iter.next())
  return outDict

def _parseSetFL(infile, finnisSinclair = False):
  """Parse text file containing setfl formatted data into a dictionary.

  @param infile Python file object containing setfl information
  @param finnisSinclair If True, then assume that file is eam.fs format and contains multiple density functions for each element
                        If False, assume only a single density function per element"""
  outdict = {}

  #Read comments
  outdict['comments'] = []
  outdict['comments'].append(infile.next()[:-1])
  outdict['comments'].append(infile.next()[:-1])
  outdict['comments'].append(infile.next()[:-1])

  #Read types line
  line = infile.next()[:-1]
  tokens= line.split()
  outdict['ntypes'] = int(tokens[0])
  outdict['types'] = tokens[1:]

  #Read numsteps and cutoffs
  line = infile.next()[:-1]
  tokens = line.split()
  nrho,drho,nr,dr,rcutoff = tokens
  outdict['nrho'] = int(nrho)
  outdict['drho'] = float(drho)
  outdict['nr'] = int(nr)
  outdict['dr'] = float(dr)
  outdict['rcutoff'] = float(rcutoff)

  #Set-up function used to read the density blocks (depending on whether we're processing a finnis sinclair file or not)
  readDensityFunction = None
  def readSetFLDensityFunction(floatit, outdict, workingdict):
    for i in xrange(outdict['nrho']):
      workingdict.setdefault('densityfunction', []).append(floatit.next())

  #... this version of the function creates 'densityfunctions' list and creates NumElements^2 lists for density functions
  def readFinnisSinclairDensityFunction(floatit, outdict, workingdict):
    nrho = outdict['nrho']
    numelements = outdict['ntypes']
    workinglist = []
    for n in xrange(numelements):
      workingdenslist = []
      for i in xrange(nrho):
        workingdenslist.append(floatit.next())
      workinglist.append(workingdenslist)
    workingdict['densityfunctions'] = workinglist

  if finnisSinclair == False:
    readDensityFunction = readSetFLDensityFunction
  else:
    readDensityFunction = readFinnisSinclairDensityFunction

  #Read embedding function and density function
  def readelementblock():
    workingdict = {}
    headerline = infile.next()[:-1]
    ielem, amass, blat, lat = headerline.split()
    workingdict['ielem'] = int(ielem)
    workingdict['amass'] = float(amass)
    workingdict['blat'] = float(blat)
    workingdict['lat'] = lat
    floatit = _floatiterator(infile)
    for i in xrange(outdict['nr']):
      workingdict.setdefault('embeddingfunction', []).append(floatit.next())
    floatit = _floatiterator(infile)
    readDensityFunction(floatit, outdict, workingdict)
    outdict.setdefault('elementblocks', []).append(workingdict)

  for i in xrange(outdict['ntypes']):
    readelementblock()

  #Read pair potentials
  def readppairblock():
    workinglist = []
    floatit = _floatiterator(infile)
    for i in xrange(outdict['nr']):
      workinglist.append(floatit.next())
    outdict.setdefault('ppairs', []).append(workinglist)

  for i in xrange(outdict['ntypes']):
    for j in xrange(i,outdict['ntypes']):
      readppairblock()
  return outdict


class RunLAMMPSEAMTableTestCase(TempfileTestCase):
  """TestCase that runs lammps if possible"""

  @needsLAMMPS
  def testFluorite_NoPair(self):
    """Test EAM Tabulation using pair potentials for a fluorite structure"""
    shutil.copyfile(
      os.path.join(_getResourceDirectory(), 'calc_energy.lmpin'),
      os.path.join(self.tempdir, 'calc_energy.lmpin'))

    shutil.copyfile(
      os.path.join(_getResourceDirectory(), 'CeO2-single_cell.lmpstruct'),
      os.path.join(self.tempdir, 'structure.lmpstruct'))

    oldpwd = os.getcwd()
    os.chdir(self.tempdir)
    try:
      with open("potentials.lmpinc", "wb") as potfile:
        print >>potfile, "pair_style eam/fs"
        print >>potfile, "pair_coeff   *    *  eam.fs O Ce"
        print >>potfile, ""
        print >>potfile, "replicate 4 4 4"


      # Now define potentials
      def erf(x):
          # save the sign of x
          sign = 1 if x >= 0 else -1
          x = abs(x)

          # constants
          a1 =  0.254829592
          a2 = -0.284496736
          a3 =  1.421413741
          a4 = -1.453152027
          a5 =  1.061405429
          p  =  0.3275911

          # A&S formula 7.1.26
          t = 1.0/(1.0 + p*x)
          y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*math.exp(-x*x)
          return sign*y # erf(-x) = -erf(x)


      def morse(gamma, r_star, D):
        def f(r):
          return D*(math.exp(-2.0*gamma*(r-r_star)) - 2.0*math.exp(-gamma*(r-r_star)))
        return f

      def zw(origfunc):
        def f(r):
          if r == 0.0:
            return 0.0
          return origfunc(r)
        return f

      def zeropot(r):
        return 0.0

      def embedCe(d):
        return -0.30766574471*(d**0.5)

      def embedO(d):
        return -0.690017059438*(d**0.5)

      def densityCe(r):
        m=8
        a=1556.80263774 #tmpa #55,100//0.00002,0.0001
        return (a/(r**m))*0.5*(1+erf(20*(r-1.5)))

      def densityO(r):
        m=8
        a=106.855913747 #tmpa #55,100//0.00002,0.0001
        return (a/(r**m))*0.5*(1+erf(20*(r-1.5)))

      functionO_Ce = zw(potentials.plus(potentials.buck(351.341192796, 0.380516580733, 0.0), morse(1.86874578949, 2.35603812582, 0.719250701502)))
      functionO_O = zw(potentials.buck(830.283447557, 0.352856254215, 3.88437209048))
      functionCe_Ce = zw(potentials.buck(18600.0, 0.26644, 0.0))

      potO_O = potentials.Potential('O', 'O', functionO_O)
      potO_Ce = potentials.Potential('O', 'Ce', functionO_Ce)
      potCe_Ce = potentials.Potential('Ce', 'Ce', functionCe_Ce)

      eampotCe = potentials.EAMPotential('Ce', 92, 238, embedCe, {'Ce': zw(densityCe), 'O' : zw(densityCe)})
      eampotO = potentials.EAMPotential('O', 8, 16, embedO, {'Ce': zw(densityO), 'O' : zw(densityO)})

      eampots = [eampotCe, eampotO]
      pairpots = [potO_O, potO_Ce, potCe_Ce]

      drho = 0.001
      nrho = int(100.0 / drho)

      cutoff = 11.0
      dr = 0.01
      nr = int(cutoff/dr)

      with open('eam.fs', 'wb') as outfile:
        potentials.writeSetFLFinnisSinclair(
            nrho, drho,
            nr, dr,
            eampots,
            [],
            outfile,
            ["Yakub potential. Zero embed and density", "",""],
            cutoff)

      runLAMMPS()
      energy = extractLAMMPSEnergy()
      expect = -1001.27446044
      self.assertAlmostEquals(expect, energy, places=5)
    finally:
      os.chdir(oldpwd)

  @needsLAMMPS
  def testWriteFuncflPair(self):
    """Test that unit conversions for pair potential tabulation by writeFuncFL() are correct"""
    shutil.copyfile(os.path.join(_getResourceDirectory(), "writefuncfl_pair.lmpstruct"), os.path.join(self.tempdir,"structure.lmpstruct"))
    shutil.copyfile(os.path.join(_getResourceDirectory(), "calc_energy.lmpin"), os.path.join(self.tempdir,"calc_energy.lmpin"))

    oldpwd = os.getcwd()
    os.chdir(self.tempdir)
    try:
      with open("potentials.lmpinc", "wb") as potfile:
        print >>potfile, "pair_style eam"
        print >>potfile, "pair_coeff 1 1 Ag.eam"

      def embed(rho):
        return 0.0

      def density(rij):
        return 0.0

      def pair_AgAg(rij):
        if rij == 0:
          return 0.0
        return (2.485883762/rij) ** 12

      pairPotentials = [ Potential('Ag', 'Ag', pair_AgAg) ]

      # Create EAMPotential
      eamPotentials = [ EAMPotential("Ag", 47, 107.8682, embed, density) ]

      nrho = 50000
      drho = 0.001

      nr = 12000
      dr = 0.001

      from atsim.potentials import writeFuncFL

      with open("Ag.eam", 'wb') as outfile:
        writeFuncFL(
          nrho, drho,
          nr, dr,
          eamPotentials,
          pairPotentials,
          title='Sutton Chen Ag',
          out= outfile)

      runLAMMPS()
      energy = extractLAMMPSEnergy()
      self.assertAlmostEquals(982.2756583, energy)
    finally:
      os.chdir(oldpwd)

  @needsLAMMPS
  def testWriteFuncflDensity(self):
    """Test writeFuncFL() writing of density function correct"""
    shutil.copyfile(os.path.join(_getResourceDirectory(), "writefuncfl_pair.lmpstruct"), os.path.join(self.tempdir,"structure.lmpstruct"))
    shutil.copyfile(os.path.join(_getResourceDirectory(), "calc_energy.lmpin"), os.path.join(self.tempdir,"calc_energy.lmpin"))
    oldpwd = os.getcwd()
    os.chdir(self.tempdir)
    try:
      with open("potentials.lmpinc", "wb") as potfile:
        print >>potfile, "pair_style eam"
        print >>potfile, "pair_coeff 1 1 Ag.eam"


      def embed(rho):
        return rho

      def density(rij):
        if rij == 0:
          return 0.0
        return (2.928323832 / rij) ** 6.0

      def pair_AgAg(rij):
        return 0.0

      pairPotentials = [ Potential('Ag', 'Ag', pair_AgAg) ]

      # Create EAMPotential
      eamPotentials = [ EAMPotential("Ag", 47, 107.8682, embed, density) ]

      nrho = 200000
      drho = 0.0005

      nr = 5000
      dr = 0.001

      from atsim.potentials import writeFuncFL

      with open("Ag.eam", 'wb') as outfile:
        writeFuncFL(
          nrho, drho,
          nr, dr,
          eamPotentials,
          pairPotentials,
          title='Sutton Chen Ag',
          out= outfile)

      runLAMMPS()
      energy = extractLAMMPSEnergy()
      self.assertAlmostEquals(83.7425918012, energy/2.0, places = 5)
    finally:
      os.chdir(oldpwd)

  @needsLAMMPS
  def testWriteFuncflEmbed(self):
    """Test that writing of embedding function by writeFuncFL() is correct"""
    shutil.copyfile(os.path.join(_getResourceDirectory(), "writefuncfl_pair.lmpstruct"), os.path.join(self.tempdir,"structure.lmpstruct"))
    shutil.copyfile(os.path.join(_getResourceDirectory(), "calc_energy.lmpin"), os.path.join(self.tempdir,"calc_energy.lmpin"))
    oldpwd = os.getcwd()
    os.chdir(self.tempdir)
    try:
      with open("potentials.lmpinc", "wb") as potfile:
        print >>potfile, "pair_style eam"
        print >>potfile, "pair_coeff 1 1 Ag.eam"


      def embed(rho):
        return -math.sqrt(rho)

      def density(rij):
        return 2.0

      def pair_AgAg(rij):
        return 0.0

      pairPotentials = [ Potential('Ag', 'Ag', pair_AgAg) ]

      # Create EAMPotential
      eamPotentials = [ EAMPotential("Ag", 47, 107.8682, embed, density) ]

      nrho = 100000
      drho = 0.001

      nr = 5000
      dr = 0.001

      from atsim.potentials import writeFuncFL

      with open("Ag.eam", 'wb') as outfile:
        writeFuncFL(
          nrho, drho,
          nr, dr,
          eamPotentials,
          pairPotentials,
          title='Sutton Chen Ag',
          out= outfile)

      runLAMMPS()
      energy = extractLAMMPSEnergy()
      self.assertAlmostEquals(-math.sqrt(2), energy/2.0)
    finally:
      os.chdir(oldpwd)


  def createSetFlFuncs(self):
    def makeFunc(a, b, r_e, c):
      def func(r):
        return (a * math.exp(-b*(r/r_e -1)))/(1+(r/r_e - c)**20.0)
      return func

    def makePairPotAA(A, gamma, r_e, kappa,
                      B, omega, lamda):
      f1 = makeFunc(A, gamma, r_e, kappa)
      f2 = makeFunc(B, omega, r_e, lamda)
      def func(r):
        return f1(r) - f2(r)
      return func

    def makePairPotAB(dens_a, phi_aa, dens_b, phi_bb):
      def func(r):
        return 0.5 * ( (dens_b(r)/dens_a(r) * phi_aa(r)) + (dens_a(r)/dens_b(r) * phi_bb(r)) )
      return func


    def makeEmbed(rho_e, rho_s, F_ni, F_i, F_e, nu):
      rho_n = 0.85*rho_e
      rho_0 = 1.15*rho_e

      def e1(rho):
        return sum([F_ni[i] * (rho/rho_n - 1)**float(i) for i in xrange(4)])

      def e2(rho):
        return sum([F_i[i] * (rho/rho_e - 1)**float(i) for i in xrange(4)])

      def e3(rho):
        return F_e * (1.0 - nu*math.log(rho/rho_s)) * (rho/rho_s)**nu

      def func(rho):
        if rho < rho_n:
          return e1(rho)
        elif rho_n <= rho < rho_0:
          return e2(rho)
        return e3(rho)
      return func

    r_eCu     = 2.556162
    f_eCu     = 1.554485
    gamma_Cu  = 8.127620
    omega_Cu  = 4.334731
    A_Cu      = 0.396620
    B_Cu      = 0.548085
    kappa_Cu  = 0.308782
    lambda_Cu = 0.756515

    rho_e_Cu  = 21.175871
    rho_s_Cu  = 21.175395
    F_ni_Cu   = [-2.170269, -0.263788, 1.088878, -0.817603]
    F_i_Cu    = [-2.19, 0.0, 0.561830, -2.100595]
    nu_Cu     = 0.310490
    F_e_Cu    = -2.186568

    r_eAl     = 2.863924
    f_eAl     = 1.403115
    gamma_Al  = 6.613165
    omega_Al  = 3.527021
    # A_Al      = 0.134873
    A_Al      = 0.314873
    B_Al      = 0.365551
    kappa_Al  = 0.379846
    lambda_Al = 0.759692

    rho_e_Al  =  20.418205
    rho_s_Al  =  23.195740
    F_ni_Al   = [-2.807602, -0.301435, 1.258562, -1.247604]
    F_i_Al    = [-2.83, 0.0, 0.622245, -2.488244]
    nu_Al     = 0.785902
    F_e_Al    = -2.824528

    dens_Cu   = makeFunc(f_eCu, omega_Cu, r_eCu, lambda_Cu )
    dens_Al   = makeFunc(f_eAl, omega_Al, r_eAl,  lambda_Al )

    pair_CuCu = makePairPotAA(A_Cu, gamma_Cu, r_eCu, kappa_Cu,
                             B_Cu, omega_Cu, lambda_Cu)

    pair_AlAl = makePairPotAA(A_Al, gamma_Al, r_eAl, kappa_Al,
                             B_Al, omega_Al, lambda_Al)

    pair_AlCu = makePairPotAB(dens_Cu, pair_CuCu, dens_Al, pair_AlAl)

    embed_Cu = makeEmbed(rho_e_Cu, rho_s_Cu, F_ni_Cu, F_i_Cu, F_e_Cu, nu_Cu)
    embed_Al = makeEmbed(rho_e_Al, rho_s_Al, F_ni_Al, F_i_Al, F_e_Al, nu_Al)


    return {
      'dens_Cu' : dens_Cu,
      'dens_Al' : dens_Al,
      'pair_CuCu' : pair_CuCu,
      'pair_AlAl' : pair_AlAl,
      'pair_AlCu' : pair_AlCu,
      'embed_Cu' : embed_Cu,
      'embed_Al' : embed_Al}


  def createSetFlPairPots(self):
    fdict = self.createSetFlFuncs()

    pair_AlAl = fdict['pair_AlAl']
    pair_CuCu = fdict['pair_CuCu']
    pair_AlCu = fdict['pair_AlCu']

    pairPotentials = [
        Potential('Al', 'Al', pair_AlAl),
        Potential('Cu', 'Cu', pair_CuCu),
        Potential('Al', 'Cu', pair_AlCu)]
    return pairPotentials

  def createSetFLEAMPots(self):
    fdict = self.createSetFlFuncs()

    eamPotentials = [
      EAMPotential("Al", 13, 26.98, fdict['embed_Al'], fdict['dens_Al']),
      EAMPotential("Cu", 29, 63.55, fdict['embed_Cu'], fdict['dens_Cu'])]

    return eamPotentials


class LAMMPSWriteEAMTableTestCase(unittest.TestCase):
  """TestCase for the lammps.writeEAMTable module"""

  def testWriteFuncflAgU3EAM(self):
    """Test that lammps.writeEAMTable can reproduce a DYNAMO funcfl EAM file from the standard LAMMPS distribution here Ag_u3.eam"""

    #Read the expected table from the .eam file
    with contextlib.closing(open(os.path.join(_getResourceDirectory(), 'Ag_u3.eam'))) as infile:
      expectTable = _parseEAMTable(infile)

    #Create the potential callables to be passed into writeEAMTable()
    embeddingFunction = potentials.TableReader( open(os.path.join(_getResourceDirectory(), 'AgU3embedding.table'), 'r'))
    effectiveChargeFunction = potentials.TableReader( open(os.path.join(_getResourceDirectory(), 'AgU3effectivecharge.table'), 'r'))

    def effectiveChargeFunction_eV(rij):
      v = effectiveChargeFunction(rij)
      v *= v
      v *= 27.2 * 0.529
      if rij != 0.0:
        v /= rij
      return v

    densityFunction = potentials.TableReader(
      open(os.path.join(_getResourceDirectory(), 'AgU3density.table'), 'r'))

    title = "Ag functions (universal 3), SM Foiles et al, PRB, 33, 7983 (1986)"
    atomicNumber = 47
    mass = 107.87
    latticeConstant = 4.09
    latticeType = 'FCC'

    nrho = 500
    drho = 5.0100200400801306e-04

    nr = 500
    dr = 1.1212121212121229e-02

    eampotlist = [EAMPotential("Ag", atomicNumber, mass, embeddingFunction, densityFunction, latticeConstant, latticeType)]
    potlist = [Potential("Ag", "Ag", effectiveChargeFunction_eV)]

    actualTable = StringIO.StringIO()
    potentials.writeFuncFL(
      nrho, drho,
      nr, dr,
      eampotlist,
      potlist,
      title = title,
      out = actualTable)
    actualTable.seek(0)
    actualTable = _parseEAMTable(actualTable)
    import testutil
    testutil.compareCollection(self,expectTable, actualTable, places = 3)


  def testWriteSetFLFromParameters(self):
    """Test that lammps.potentials.writeSetFL() can generate the Al_zhou.eam.alloy file from LAMMPS distribution"""

    with contextlib.closing(open( os.path.join(_getResourceDirectory(), 'Al_zhou.eam.alloy'), 'r')) as infile:
      expectEAMTable = _parseSetFL(infile)

    comments = [
      '#-> LAMMPS Potential File in DYNAMO 86 setfl Format <-#',
      '# Zhou Al Acta mater(2001)49:4005',
      '# Implemented by G. Ziegenhain (2007) gerolf@ziegenhain.com']

    #Density function
    def densityFunction(r):
      beta = 3.702623
      r_e  = 2.886166
      lamb = 0.790264
      f_e  = 1.392302
      num = f_e * math.exp( -beta * ( (r/r_e) - 1.0) )
      den = 1.0 + (( (r/r_e) - lamb )**20)
      return num/den

    #Embedding function
    def subembed(i, F_ni, rho_n, rho):
      return F_ni * ((rho/rho_n) - 1.0)**float(i)

    def embeddingFunction(rho):
      rho_e = 20.226537
      rho_n = 0.85 * rho_e
      rho_o = 1.15 * rho_e
      nu = 0.779208
      if rho < rho_n:
        return subembed(0, -2.806783, rho_n, rho) + \
               subembed(1, -0.276173, rho_n, rho) + \
               subembed(2,  0.893409, rho_n, rho) + \
               subembed(3, -1.637201, rho_n, rho)
      elif rho_n <= rho < rho_o:
        return subembed(0, -2.806783, rho_e, rho) + \
               subembed(1, -0.276173, rho_e, rho) + \
               subembed(2,  0.893409, rho_e, rho) + \
               subembed(3, -1.637201, rho_e, rho)
      else:
        return -2.829437 * ( 1.0 - math.log( (rho/rho_e)**nu)) * (rho/rho_e)**nu

    #Pair potential
    def ppotfunc(r):
      A = 0.251519
      alpha = 6.94219
      r_e = 2.886166
      kappa = 0.395132
      B =   0.313394
      beta = 3.702623
      lamb = 0.790264
      return ((A * math.exp( -alpha* ( (r/r_e) - 1.0) ))/ (1.0 + ( (r/r_e - kappa)**20))) - \
             ((B * math.exp( -beta * ( (r/r_e) - 1.0) ))/ (1.0 + ( (r/r_e - lamb)**20)))
    alpp = potentials.Potential('Al', 'Al', ppotfunc)

    nrho =  10001
    drho =  0.00559521603477821424
    nr   =  10001
    dr   =  0.00101014898510148996
    cutoff = 10.10250000000000092371
    eampots = [
      potentials.EAMPotential('Al', 1, 26.982,  embeddingFunction, densityFunction, latticeConstant = 4.05,  latticeType = 'FCC')]

    pairpots = [ alpp ]

    actualEAMTable = StringIO.StringIO()
    potentials.writeSetFL(
      nrho, drho,
      nr, dr,
      eampots,
      pairpots,
      comments = comments,
      out = actualEAMTable,
      cutoff = cutoff )
    actualEAMTable.seek(0)
    actualEAMTable = _parseSetFL(actualEAMTable)

    import testutil
    testutil.compareCollection(self,expectEAMTable, actualEAMTable, places = 3, percenttolerance = 7.0)

  def testWriteSetFL(self):
    """Test creation of DYNAMO setfl formatted file for use with lammps pair_style eam/alloy"""
    with contextlib.closing(open( os.path.join(_getResourceDirectory(), 'AlCu.eam.alloy'), 'r')) as infile:
      expectEAMTable = _parseSetFL(infile)

    comments = [
    '##Al-Cu EAM potentials from J. Cai and Y.Y. Ye##',
    '##Phys. Rev. B 54, 8398-8410 (1996).############',
    '################################################']

    alEmbedFunc = TableReader(_openResource('setfl_AlEmbed.table'))
    alElectronDensity = TableReader(_openResource('setfl_AlDensity.table'))
    cuEmbedFunc = TableReader(_openResource('setfl_CuEmbed.table'))
    cuElectronDensity = TableReader(_openResource('setfl_CuDensity.table'))

    pairpot_Al_Al = TableReader(_openResource('setfl_AlAlPair.table'))
    pairpot_Cu_Al = TableReader(_openResource('setfl_CuAlPair.table'))
    pairpot_Cu_Cu = TableReader(_openResource('setfl_CuCuPair.table'))

    eampots = [
      potentials.EAMPotential('Al', 13, 26.982,  alEmbedFunc, alElectronDensity, latticeConstant = 4.05,  latticeType = 'FCC'),
      potentials.EAMPotential('Cu', 29, 63.546, cuEmbedFunc, cuElectronDensity, latticeConstant = 3.615, latticeType = 'FCC')]

    pairpots = [
      Potential('Al', 'Al', pairpot_Al_Al),
      Potential('Cu', 'Al', pairpot_Cu_Al),
      Potential('Cu', 'Cu', pairpot_Cu_Cu) ]

    nrho = 1000
    drho = 1.0396723078e+00

    nr = 3000
    dr = 2.2282427476e-03

    cutoff = 6.6825000000e+00

    actualEAMTable = StringIO.StringIO()
    potentials.writeSetFL(
      nrho, drho,
      nr, dr,
      eampots,
      pairpots,
      comments = comments,
      out = actualEAMTable,
      cutoff = cutoff)
    actualEAMTable.seek(0)
    actualEAMTable = _parseSetFL(actualEAMTable)

    import testutil
    testutil.compareCollection(self,expectEAMTable, actualEAMTable, places = 3)

  def testWriteSetFLFinnisSinclair(self):
    """Test that lammps.potentials.writeSetFLFinnisSinclair() (suitable for use with pair_style eam/fs) can re-create AlFe_mm.eam.fs file from lammps distribution"""

    #Open the expected output
    with contextlib.closing(_openResource('AlFe_mm.eam.fs')) as infile:
      expectEAMTable = _parseSetFL(infile, finnisSinclair = True)

    #Set-up tabulation of the actual EAM table from parameters found in
    #M.I. Mendelev,  D.J. Srolovitz,  G.J. Ackland and  S. Han, J. Mater. Res. 20, 208-218 (2005).

    comments = [
      "Sourse: M.I. Mendelev,  D.J. Srolovitz,  G.J. Ackland and  S. Han, J. Mater. Res. 20, 208-218 (2005).",
      "Contact information: mendelev@ameslab.gov",
      "Sunday, Jun 10, 2007  The potential was taken from Al3Fe_D03 (in C:\SIMULATION.MD\Al-Fe\T=0)" ]

    #Embedding functions
    def alEmbedFunction(rho):
      if rho == 0.0: return 0.0
      return -math.sqrt(rho)+0.000093283590195398*rho**2-0.0023491751192724*rho*math.log(rho)

    def feEmbedFunction(rho):
      return -math.sqrt(rho) -  0.00067314115586063*rho**2 +  0.000000076514905604792*rho**4

    #Density functions
    def alAlDensFunction(r):
      funcs = [
        (2.5, lambda r: 0.00019850823042883 * (2.5 - r)**4 ),
        (2.6, lambda r: 0.10046665347629    * (2.6 - r)**4 ),
        (2.7, lambda r: 1.0054338881951E-01 * (2.7 - r)**4 ),
        (2.8, lambda r: 0.099104582963213   * (2.8 - r)**4 ),
        (3.0, lambda r: 0.090086286376778   * (3.0 - r)**4 ),
        (3.4, lambda r: 0.0073022698419468  * (3.4 - r)**4 ),
        (4.2, lambda r: 0.014583614223199   * (4.2 - r)**4 ),
        (4.8, lambda r: -0.0010327381407070 * (4.8 - r)**4 ),
        (5.6, lambda r: 0.0073219994475288  * (5.6 - r)**4 ),
        (6.5, lambda r: 0.0095726042919017  * (6.5 - r)**4 )]
      vals = [ func(r) for cutoff, func in funcs if r <= cutoff ]
      return sum(vals)

    def feFeDensFunction(r):
      funcs = [
        (2.4, lambda r: 11.686859407970*(2.4    - r)**3),
        (3.2, lambda r: -0.014710740098830*(3.2 - r)**3),
        (4.2, lambda r: 0.47193527075943*(4.2   - r)**3)]
      vals = [ func(r) for cutoff, func in funcs if r <= cutoff ]
      return sum(vals)

    def feAlDensFunction(r):
      funcs = [
        (2.4, lambda r: 0.010015421408039   * (2.4 - r)**4),
        (2.5, lambda r: 0.0098878643929526  * (2.5 - r)**4),
        (2.6, lambda r: 0.0098070326434207  * (2.6 - r)**4),
        (2.8, lambda r: 0.0084594444746494 * (2.8 - r)**4),
        (3.1, lambda r: 0.0038057610928282  * (3.1 - r)**4),
        (5.0, lambda r: -0.0014091094540309 * (5.0 - r)**4),
        (6.2, lambda r: 0.0074410802804324  * (6.2 - r)**4)]
      vals = [ func(r) for cutoff, func in funcs if r <= cutoff ]
      return sum(vals)

    def zerowrap(wrapped):
      def f(r):
        if r == 0.0:
          return 0.0
        return wrapped(r)
      return f

    #Pair potentials
    def ppfuncAlAl(r):


      funcs = [
        ( (0.0, 1.6) , zerowrap(lambda r: (2433.5591473227/r)  * (0.1818 * math.exp(-22.713109144730* r) + 0.5099 * math.exp(-6.6883008584622* r) + 0.2802* math.exp(-2.8597223982536*r) + 0.02817* math.exp(-1.4309258761180*r)) )),
        ( (1.6, 2.25), lambda r: math.exp(6.0801330531321 - 2.3092752322555*r + 0.042696494305190*r**2 - 0.07952189194038*r**3) ),
        ( (2.25, 3.2), lambda r:  17.222548257633*(3.2 - r)**4 - 13.838795389103*(3.2 - r)**5  + 26.724085544227*(3.2 - r)**6 - 4.8730831082596*(3.2 - r)**7 + 0.26111775221382*(3.2 - r)**8),
        ( (2.25, 4.8), lambda r: -1.8864362756631*(4.8 - r)**4 + 2.4323070821980*(4.8 - r)**5 - 4.0022263154653*(4.8 - r)**6 + 1.3937173764119*(4.8 - r)**7 - 0.31993486318965*(4.8 - r)**8),
        ( (2.25, 6.5), lambda r: 0.30601966016455*(6.5 - r)**4 - 0.63945082587403*(6.5 - r)**5 + 0.54057725028875*(6.5 - r)**6 - 0.21210673993915*(6.5 - r)**7 + 0.032014318882870*(6.5 - r)**8) ]

      vals = [ func(r) for ((lowcut, highcut), func) in funcs if lowcut <= r < highcut ]
      return sum(vals)

    def ppfuncAlFe(r):
      funcs = [
        ( (0.0, 1.2), zerowrap(lambda r: (4867.1182946454/r) * (0.1818 * math.exp(-25.834107666296*r) + 0.5099*math.exp(-7.6073373918597*r) + 0.2802*math.exp(-3.2526756183596*r) + 0.02817*math.exp(-1.6275487829767*r))) ),
        ( (1.2, 2.2), lambda r: math.exp(6.6167846784367-1.5208197629514*r - 0.73055022396300*r**2 - 0.038792724942647*r**3) ),
        ( (2.2, 3.2), lambda r: -4.1487019439249*(3.2 - r)**4 + 5.6697481153271*(3.2 - r)**5 - 1.7835153896441*(3.2 - r)**6 - 3.3886912738827*(3.2 - r)**7 + 1.9720627768230*(3.2 - r)**8),
        ( (2.2, 6.2), lambda r: 0.094200713038410*(6.2 - r)**4 -0.16163849208165*(6.2 - r)**5 + 0.10154590006100*(6.2 - r)**6 -0.027624717063181*(6.2 - r)**7 + 0.0027505576632627*(6.2 - r)**8) ]
      vals = [ func(r) for ((lowcut, highcut), func) in funcs if lowcut <= r < highcut ]
      return sum(vals)

    def ppfuncFeFe(r):
      funcs = [
      ( (0.0, 1.0), zerowrap(lambda r: (9734.2365892908/r)*(0.1818*math.exp(-28.616724320005*r) + 0.5099*math.exp(-8.4267310396064*r) + 0.2802*math.exp(-3.6030244464156*r) + 0.02817*math.exp(-1.8028536321603*r))) ),
      ( (1.0, 2.05), lambda r: math.exp(7.4122709384068 - 0.64180690713367*r - 2.6043547961722*r**2 + 0.62625393931230*r**3) ),
      ( (2.05, 2.2), lambda r: -27.444805994228*(2.2 - r)**3 ),
      ( (2.05, 2.3), lambda r: 15.738054058489*(2.3 - r)**3),
      ( (2.05, 2.4), lambda r: 2.2077118733936*(2.4 - r)**3),
      ( (2.05, 2.5), lambda r: -2.4989799053251*(2.5 - r)**3),
      ( (2.05, 2.6), lambda r: 4.2099676494795*(2.6 - r)**3),
      ( (2.05, 2.7), lambda r: -0.77361294129713*(2.7 - r)**3),
      ( (2.05, 2.8), lambda r: 0.80656414937789*(2.8 - r)**3),
      ( (2.05, 3.0), lambda r: -2.3194358924605*(3.0 - r)**3),
      ( (2.05, 3.3), lambda r: 2.6577406128280*(3.3 - r)**3),
      ( (2.05, 3.7), lambda r: -1.0260416933564*(3.7 - r)**3),
      ( (2.05,4.2), lambda r: 0.35018615891957*(4.2 - r)**3),
      ( (2.05, 4.7), lambda r: -0.058531821042271*(4.7 - r)**3),
      ( (2.05, 5.3), lambda r: -0.0030458824556234*(5.3 - r)**3) ]
      vals = [ func(r) for ((lowcut, highcut), func) in funcs if lowcut <= r < highcut ]
      return sum(vals)

    class ErrorPotential(potentials.Potential):

      def energy(self, r):
        #There seems to be an error in the way that the lammps potential is tabulated (there is a step function at the start of each pair potential"
        #the following if statement has been introduced to reproduce this error and allow the test to pass
        if r <= 0.500:
          if r != 0.0:
            return 1e12/r
          else:
            return 1e12
        return potentials.Potential.energy(self, r)

    pairpots = [ ErrorPotential('Al', 'Al', ppfuncAlAl),
                 ErrorPotential('Al', 'Fe', ppfuncAlFe),
                 ErrorPotential('Fe', 'Fe', ppfuncFeFe)]

    #Other parameters
    nrho = 10000
    drho = 3.00000000000000E-0002
    nr   = 10000
    dr   = 6.50000000000000E-0004
    cutoff = 6.5

    #Assemble the EAMPotential objects
    eampots = [
      #Al
      potentials.EAMPotential('Al', 13, 26.98154, alEmbedFunction,
        { 'Al' : alAlDensFunction,
          'Fe' : feAlDensFunction },
        latticeConstant = 4.04527,
        latticeType = 'fcc'),
      #Fe
      potentials.EAMPotential('Fe', 26, 55.845, feEmbedFunction,
        { 'Al': feAlDensFunction,
          'Fe' : feFeDensFunction},
        latticeConstant = 2.855312,
        latticeType = 'bcc') ]

    #Now actually generate the actual tabulated potential

    actualEAMTable = StringIO.StringIO()
    potentials.writeSetFLFinnisSinclair(
      nrho, drho,
      nr, dr,
      eampots,
      pairpots,
      comments = comments,
      out = actualEAMTable,
      cutoff = cutoff)
    actualEAMTable.seek(0)
    actualEAMTable = _parseSetFL(actualEAMTable, finnisSinclair = True)

    import testutil
    testutil.compareCollection(self,expectEAMTable, actualEAMTable, places = 3, percenttolerance = 7.0)

  def testWriteSetFLFinnisSinclair_PairPotentialOrder(self):
    """Check that pair potentials are written in the correct order"""
    nrho = 5
    drho = 0.1
    nr = 5
    dr = 0.1
    cutoff = 0.5
    import collections
    def defaultdens():
      return lambda x: 0.0

    defaultdict = collections.defaultdict(defaultdens)

    eampot1 =potentials.EAMPotential('A', 1, 1.0, lambda x: 0.1, defaultdict)
    eampot2 =potentials.EAMPotential('B', 2, 2.0, lambda x: 0.2, defaultdict)
    eampot3 =potentials.EAMPotential('C', 3, 3.0, lambda x: 0.3, defaultdict)

    pairpot_aa =potentials.Potential('A', 'A', lambda r: 1.0)
    pairpot_bb =potentials.Potential('B', 'B', lambda r: 2.0)
    pairpot_cc =potentials.Potential('C', 'C', lambda r: 3.0)
    pairpot_ba =potentials.Potential('B', 'A', lambda r: 5.0)
    pairpot_ac =potentials.Potential('A', 'C', lambda r: 6.0)
    pairpot_bc =potentials.Potential('B', 'C', lambda r: 7.0)

    # Define two species
    sio = StringIO.StringIO()
    potentials.writeSetFLFinnisSinclair(
      nrho, drho,
      nr, dr,
      [eampot2, eampot1],
      [pairpot_bb, pairpot_aa, pairpot_ba],
      comments = ['Comment1', 'Comment2', 'Comment3'],
      out = sio,
      cutoff = cutoff)

    rvals = [0.0, 0.1, 0.2, 0.3, 0.4]

    expect = [ r * pairpot_bb.energy(r) for r in rvals]
    expect.extend([ r * pairpot_ba.energy(r) for r in rvals])
    expect.extend([ r * pairpot_aa.energy(r) for r in rvals])

    sio.seek(0)
    actual = sio.readlines()
    actual = actual[37:]
    actual = [float(v) for v in actual]

    import testutil
    testutil.compareCollection(self,expect, actual)

    # Try a ternary system
    # Define two species
    sio = StringIO.StringIO()
    potentials.writeSetFLFinnisSinclair(
      nrho, drho,
      nr, dr,
      [eampot2, eampot1, eampot3],
      [pairpot_bb, pairpot_aa, pairpot_ba, pairpot_cc, pairpot_ac, pairpot_bc],
      sio,
      ['Comment1', 'Comment2', 'Comment3'],
      cutoff)

    rvals = [0.0, 0.1, 0.2, 0.3, 0.4]

    expect = [ r * pairpot_bb.energy(r) for r in rvals]
    expect.extend([ r * pairpot_ba.energy(r) for r in rvals])
    expect.extend([ r * pairpot_aa.energy(r) for r in rvals])
    expect.extend([ r * pairpot_bc.energy(r) for r in rvals])
    expect.extend([ r * pairpot_ac.energy(r) for r in rvals])
    expect.extend([ r * pairpot_cc.energy(r) for r in rvals])

    sio.seek(0)
    actual = sio.readlines()
    actual = actual[68:]
    actual = [float(v) for v in actual]

    testutil.compareCollection(self,expect, actual)
