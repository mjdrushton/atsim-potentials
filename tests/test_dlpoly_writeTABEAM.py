import unittest

from atsim import potentials
from atsim.potentials import _dlpoly_writeTABEAM

import math
import shutil
import os

from _tempfiletestcase import TempfileTestCase
from _rundlpoly import needsDLPOLY, extractDLPOLYEnergy, runDLPoly


def _getResourceDirectory():
  """Returns path to resources used by this test module (currently assumed to be sub-directory
  of test module called resources)"""
  return os.path.join(os.path.dirname(__file__), 'dl_poly_resources')

def _densityFunction(A, n, B, r0, m):
  def func(r):
    return A * r**n * math.exp(-B * (r - r0)**m)
  return func

def _stillingerWeber(A, B, rho, rmax):
  def func(r):
    if r >= rmax:
      return 0.0
    return A * math.exp(rho/(r-rmax)) * ((B/r**4) - 1.0)
  return func

def _zerowrap(func):
  def f(r):
    if r == 0.0:
      return 0.0
    return func(r)
  return f

def embeddingFunction(rho):
  return -math.sqrt(rho)

class DLPOLYWriteTABEAMTestCase(unittest.TestCase):
  """Tests for the _writeTABEAM module"""

  def _countSections(self, infile, sectionName):
    infile.seek(0)

    count = 0

    for line in infile:
      if line.startswith(sectionName):
        count +=1
    return count


  def testWriteEAMFile(self):
    """Test _writeTABEAM.writeTABEAM() function"""

    import StringIO

    #Define potentials
    ppbucku  = potentials.Potential('U', 'U', _zerowrap(potentials.buck(668.546808, 0.408333, 0.0)))
    ppbuckuc = potentials.Potential('U', 'C', _zerowrap(potentials.buck(30.885011, 0.814952, 0.0)))
    ppswcc   = potentials.Potential('C', 'C', _zerowrap(_stillingerWeber(1.1, 2.078, 1.368, 2.5257)))

    densU = _densityFunction( 1.301894, 2.0, 0.668659, 1.862363, 2.0)
    densC = _densityFunction(33.446287, 2.0, 1.318078, 1.512686, 2.0)

    UEAMPot = potentials.EAMPotential('U', 92, 19.05, embeddingFunction, densU )
    CEAMPot = potentials.EAMPotential('C',  6,  2.267,embeddingFunction, densC )

    maxrho = 410.0
    cutoff = 12.0
    drho = dr = 0.01
    nrho = int(maxrho/drho)
    nr = int(cutoff/dr)

    outfile = StringIO.StringIO()

    _dlpoly_writeTABEAM.writeTABEAM(
      nrho, drho,
      nr, dr,
      [UEAMPot, CEAMPot],
      [ppbucku, ppbuckuc, ppswcc],
      outfile,
      "Chartier/Van Brutzel Potentials")

    #Check structure of the file
    self.assertEquals(2, self._countSections(outfile, 'embe'))
    self.assertEquals(2, self._countSections(outfile, 'dens'))
    self.assertEquals(3, self._countSections(outfile, 'pair'))


  def testTabulateFunction(self):
    import StringIO

    testfile = StringIO.StringIO()

    def testFunction(r):
      return 2.0*r + 3.0

    from atsim.potentials._dlpoly_writeTABEAM import _tabulateFunction
    _tabulateFunction(testfile, testFunction, 5, 1.0)
    testfile.flush()
    testfile.seek(0)

    actualvalues = []
    for line in testfile:
      line = line[:-1]
      for token in line.split():
        actualvalues.append(float(token))

    expectvalues = [3.0, 5.0, 7.0, 9.0, 11.0]

    self.assertEquals(expectvalues, actualvalues)


class DLPOLYWriteTABEAMTestCase_RunDLPoly(TempfileTestCase):
  """Class that makes EAM tabulation for DLPOLY then runs code to check that we get expected energy"""

  def setUp(self):
    super(DLPOLYWriteTABEAMTestCase_RunDLPoly, self).setUp()
    shutil.copyfile(os.path.join(_getResourceDirectory(), 'CONTROL_CeO2'), os.path.join(self.tempdir, 'CONTROL'))
    shutil.copyfile(os.path.join(_getResourceDirectory(), 'CONFIG_CeO2'), os.path.join(self.tempdir, 'CONFIG'))
    shutil.copyfile(os.path.join(_getResourceDirectory(), 'FIELD_CeO2'), os.path.join(self.tempdir, 'FIELD'))

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
      return (a/(r**m))#*0.5*(1+erf(20*(r-1.5)))

    def densityO(r):
      m=8
      a=106.855913747 #tmpa #55,100//0.00002,0.0001
      return (a/(r**m))#*0.5*(1+erf(20*(r-1.5)))

    functionO_Ce = zw(potentials.plus(potentials.buck(351.341192796, 0.380516580733, 0.0), morse(1.86874578949, 2.35603812582, 0.719250701502)))
    functionO_O = zw(potentials.buck(830.283447557, 0.352856254215, 3.88437209048))
    functionCe_Ce = zw(potentials.buck(18600.0, 0.26644, 0.0))

    potO_O = potentials.Potential('O', 'O', functionO_O)
    potO_Ce = potentials.Potential('O', 'Ce', functionO_Ce)
    potCe_Ce = potentials.Potential('Ce', 'Ce', functionCe_Ce)

    eampotCe = potentials.EAMPotential('Ce', 92, 238, embedCe, zw(densityCe))
    eampotO = potentials.EAMPotential('O', 8, 16, embedO, zw(densityO))

    self.eampots = [eampotCe, eampotO]
    self.pairpots = [potO_O, potO_Ce, potCe_Ce]


  @needsDLPOLY
  def testFluorite_NoPair(self):
    # Set-up the temporary directory
    oldpwd = os.getcwd()
    os.chdir(self.tempdir)
    try:

      # Tabulate potentials
      drho = 0.001
      nrho = int(200.0 / drho)

      cutoff = 5.001
      dr = 0.001
      nr = int(cutoff/dr)

      with open('TABEAM', 'wb') as outfile:
        _dlpoly_writeTABEAM.writeTABEAM(
            nrho, drho,
            nr, dr,
            self.eampots,
            [],
            outfile,
            "Title")

      runDLPoly()
      # import pdb;pdb.set_trace()
      engcfg = extractDLPOLYEnergy()
      expect = -998.811
      self.assertAlmostEquals(expect, engcfg, places = 3)
    finally:
      os.chdir(oldpwd)

  @needsDLPOLY
  def testFluorite_WithPair(self):
    # Set-up the temporary directory
    oldpwd = os.getcwd()
    os.chdir(self.tempdir)
    try:

      # Tabulate potentials
      drho = 0.001
      nrho = int(50.0 / drho) + 1

      cutoff = 5.0
      dr = 0.001
      nr = int(cutoff/dr) +1

      with open('TABEAM', 'wb') as outfile:
        _dlpoly_writeTABEAM.writeTABEAM(
            nrho, drho,
            nr, dr,
            self.eampots,
            self.pairpots,
            outfile,
            "Title")

      runDLPoly()
      engcfg = extractDLPOLYEnergy()
      expect =-585.4435
      self.assertAlmostEquals(expect, engcfg, places = 3)
    finally:
      os.chdir(oldpwd)


class DLPOLYWriteTABEAMFinnisSinclair(TempfileTestCase):
  """Tests for the EEAM variant of the DL_POLY TABEAM file"""

  def setUp(self):
    super(DLPOLYWriteTABEAMFinnisSinclair, self).setUp()

    shutil.copyfile(os.path.join(_getResourceDirectory(),"CONFIG_triplet_eeam"), os.path.join(self.tempdir,"CONFIG"))
    shutil.copyfile(os.path.join(_getResourceDirectory(),"CONTROL_triplet_eeam"), os.path.join(self.tempdir,"CONTROL"))
    shutil.copyfile(os.path.join(_getResourceDirectory(),"FIELD_triplet_eeam"), os.path.join(self.tempdir,"FIELD"))

  @needsDLPOLY
  def testPairPotentials(self):
    """Test tabulation of pair potentials"""

    def embed(rho):
      return 0.0

    def density(r):
      return 1.0

    def ppAA(r):
      return 1.0

    def ppAB(r):
      return 2.0 * r

    def ppBB(r):
      return 9.0 * r

    Potential = potentials.Potential
    EAMPotential = potentials.EAMPotential
    pairPotentials = [
      Potential("Ar", "Ar", ppAA),
      Potential("Ar", "B", ppAB),
      Potential("B", "B", ppBB)
    ]

    eamPotentials = [
      EAMPotential("Ar", 18, 39.948, embed,
        {"Ar" : density, "B" : density}),
      EAMPotential("B", 5, 10.811, embed,
        {"Ar" : density, "B" : density})
    ]

    oldpwd = os.getcwd()
    try:
      os.chdir(self.tempdir)

      with open("TABEAM", "wb") as outfile:
        potentials.writeTABEAMFinnisSinclair(1000, 0.1, 1000, 0.01, eamPotentials, pairPotentials, outfile, "")

      runDLPoly()
      # import pdb;pdb.set_trace()
      energy = extractDLPOLYEnergy()
      self.assertAlmostEquals(65.31152, energy, places=4)

    finally:
      os.chdir(oldpwd)


  @needsDLPOLY
  def testEmbeddingFunction(self):
    """Test tabulation of embedding functions"""

    def embedA(rho):
      return 3.1 * rho

    def embedB(rho):
      return 2.13 * rho

    def density(r):
      return 1.0

    def zero(r):
      return 0.0

    Potential = potentials.Potential
    EAMPotential = potentials.EAMPotential
    pairPotentials = [
      Potential("Ar", "Ar", zero),
      Potential("Ar", "B", zero),
      Potential("B", "B", zero)
    ]

    eamPotentials = [
      EAMPotential("Ar", 18, 39.948, embedA,
       {"Ar" : density, "B" : density}),
      EAMPotential("B", 5, 10.811, embedB,
       {"Ar" : density, "B" : density})
    ]

    oldpwd = os.getcwd()
    try:
      os.chdir(self.tempdir)

      with open("TABEAM", "wb") as outfile:
        potentials.writeTABEAMFinnisSinclair(1000, 0.1, 1000, 0.01, eamPotentials, pairPotentials, outfile)

      runDLPoly()
      energy = extractDLPOLYEnergy()
      self.assertAlmostEquals(14.72, energy, places=4)

    finally:
      os.chdir(oldpwd)



  @needsDLPOLY
  def testDensityFunctions(self):
    """Test tabulation of  density functions"""

    def embed(rho):
      return rho

    def densityAA(r):
      return 1.0

    def densityBA(r):
      return 0.567*r

    def densityAB(r):
      return 0.11*r

    def densityBB(r):
      return 0.98*r


    def zero(r):
      return 0.0

    Potential = potentials.Potential
    EAMPotential = potentials.EAMPotential
    pairPotentials = [
      Potential("Ar", "Ar", zero),
      Potential("Ar", "B", zero),
      Potential("B", "B", zero)
    ]

    eamPotentials = [
      EAMPotential("Ar", 18, 39.948, embed,
       {"Ar" : densityAA, "B" : densityAB}),
      EAMPotential("B", 5, 10.811, embed,
       {"Ar" : densityBA, "B" : densityBB})
    ]

    oldpwd = os.getcwd()
    try:
      os.chdir(self.tempdir)

      with open("TABEAM", "wb") as outfile:
        potentials.writeTABEAMFinnisSinclair(1000, 0.1, 1000, 0.01, eamPotentials, pairPotentials, outfile)

      runDLPoly()
      energy = extractDLPOLYEnergy()

      expect = 0.11*5 + 0.11*2.5
      expect += 0.567*5.0 + 0.98*5.590169
      expect += 0.567*2.5 + 0.98*5.590169

      self.assertAlmostEquals(expect, energy, places=4)

    finally:
      os.chdir(oldpwd)



