"""Tests to make sure tests given in the documentation work"""

import os
import unittest
import imp
import shutil

from atsim.potentials import Potential, EAMPotential
import atsim.potentials

from _tempfiletestcase import TempfileTestCase
from _runlammps import needsLAMMPS, extractLAMMPSEnergy, runLAMMPS
from _rundlpoly import needsDLPOLY, extractDLPOLYEnergy, runDLPoly

def _getDocsDirectory():
  """Returns absolute path to docs/ directory"""
  currdir = os.path.join(os.path.dirname(__file__))
  docsdir = os.path.join("docs", "potentials")
  return os.path.abspath(docsdir)

def _getLAMMPSResourceDirectory():
  return os.path.join(os.path.dirname(__file__), 'lammps_resources')

def _getDLPolyResourceDirectory():
  return os.path.join(os.path.dirname(__file__), 'dl_poly_resources')


def _loadModule(scriptname):
  name = os.path.basename(scriptname)
  name = os.path.splitext(name)[0]
  with open(scriptname) as infile:
    mod = imp.load_module(name,infile, scriptname, ('.py', 'U', 1))
  return mod

#basak_tabulate.py
class basak_tabulateTestCase(TempfileTestCase):
  """Test docs/potentials/basak_tabulate.py"""
  test_name = os.path.join(_getDocsDirectory(), os.pardir, "basak_tabulate.py")

  @needsDLPOLY
  def testExample(self):
    """Test example basak_tabulate.py"""
    from atsim.potentials import potentialforms
    import itertools
    exampleModule = _loadModule(self.test_name)

    dldir = _getDLPolyResourceDirectory()

    oldpwd = os.getcwd()
    try:
      os.chdir(self.tempdir)


      shutil.copyfile(os.path.join(dldir, "CONTROL_pair"), "CONTROL")
      exampleModule.main()

      def dotest(d, testfunc):

        # exampleModule.main()

        pobjs = exampleModule.makePotentialObjects()
        species = [d["speciesA"], d["speciesB"]]
        spairs = [tuple(sorted(p)) for p in itertools.combinations_with_replacement(species, 2)]
        species = set(spairs)

        pobjs = [p for p in pobjs if tuple(sorted([p.speciesA, p.speciesB])) in species]

        with open('TABLE', 'wb') as outfile:
            atsim.potentials.writePotentials(
               'DL_POLY',
               pobjs,
               6.5, 6500,
               out = outfile)

        d = dict(d)
        for i in range(4):
          r = float(i)
          r += 1.0
          expect = testfunc(r)

          d['Ax'] = r

          with open(os.path.join(dldir, "CONFIG_pair.in")) as infile:
            with open("CONFIG", "wb") as outfile:
              outfile.write(infile.read() % d)


          import StringIO
          sio = StringIO.StringIO()
          print >>sio, "vdw %d" % len(pobjs)
          for p in pobjs:
            print >>sio, "%s %s tab" % (p.speciesA, p.speciesB)
          d["potDef"] = sio.getvalue()

          with open(os.path.join(dldir, "FIELD_pair.in")) as infile:
            with open("FIELD", "wb") as outfile:
              outfile.write(infile.read() % d)

          runDLPoly()
          dlenergy = extractDLPOLYEnergy()
          self.assertAlmostEqual(expect, dlenergy, places = 4)
          os.remove("STATIS")
          os.remove("CONFIG")
          os.remove("FIELD")

      d = dict(speciesA = "O", speciesB="O",
               Ax = 0.0, Ay = 0.0, Az = 0.0,
               Bx = 0.0, By = 0.0, Bz = 0.0)
      testfunc = potentialforms.buck(1633.00510, 0.327022, 3.948790)
      dotest(d, testfunc)

      d = dict(speciesA = "U", speciesB="U",
               Ax = 0.0, Ay = 0.0, Az = 0.0,
               Bx = 0.0, By = 0.0, Bz = 0.0)
      testfunc = potentialforms.buck(294.640000, 0.327022, 0.0)
      dotest(d, testfunc)

      d = dict(speciesA = "O", speciesB="U",
               Ax = 0.0, Ay = 0.0, Az = 0.0,
               Bx = 0.0, By = 0.0, Bz = 0.0)
      buck_OU = potentialforms.buck(693.648700, 0.327022, 0.0)
      morse_OU = potentialforms.morse( 1.6500, 2.36900, 0.577190)
      testfunc = atsim.potentials.plus(buck_OU, morse_OU)

      dotest(d, testfunc)

    finally:
      os.chdir(oldpwd)

  def testLAMMPSExample(self):
    """Test doumentation example Quick-Start: LAMMPS."""

    test_name = os.path.join(_getDocsDirectory(), os.pardir, "basak_tabulate_lammps.py")

    oldpwd = os.getcwd()
    try:
      os.chdir(self.tempdir)
      exampleModule = _loadModule(test_name)
      exampleModule.main()
    finally:
      os.chdir(oldpwd)



#eam_tabulate_example1.py
class eam_tabulate_example1TestCase(TempfileTestCase):
  """Test docs/potentials/eam_tabulate_example1.py"""
  test_name = os.path.join(_getDocsDirectory(), "eam_tabulate_example1.py")

  @needsLAMMPS
  def testExample(self):
    """Test example eam_tabulate_example1.py"""
    exampleModule = _loadModule(self.test_name)

    shutil.copyfile(os.path.join(_getLAMMPSResourceDirectory(), "doc_example1_Ag_fcc.lmpstruct"), os.path.join(self.tempdir,"structure.lmpstruct"))
    shutil.copyfile(os.path.join(_getLAMMPSResourceDirectory(), "calc_energy.lmpin"), os.path.join(self.tempdir,"calc_energy.lmpin"))

    oldpwd = os.getcwd()
    os.chdir(self.tempdir)
    try:
      with open("potentials.lmpinc", "wb") as potfile:
        print >>potfile, "pair_style eam"
        print >>potfile, "pair_coeff 1 1 Ag.eam"

      # Run the main method
      exampleModule.main()

      runLAMMPS()
      energy = extractLAMMPSEnergy()

      self.assertAlmostEquals(-8.23982879, energy, places = 5)
    finally:
      os.chdir(oldpwd)


class eam_tabulate_example2TestCase(TempfileTestCase):
  """Test docs/potentials/eam_tabulate_example2a.py and docs/potentials/eam_tabulate_example2b.py"""
  test_nameA = os.path.join(_getDocsDirectory(), "eam_tabulate_example2a.py")
  test_nameB= os.path.join(_getDocsDirectory(), "eam_tabulate_example2b.py")

  @needsLAMMPS
  def testExampleA_Pair(self):
    """Test pair-potentials correctly defined for EAM tabulation documentation example 2a"""
    exampleModule = _loadModule(self.test_nameA)

    shutil.copyfile(os.path.join(_getLAMMPSResourceDirectory(), "setfl_pair.lmpstruct"), os.path.join(self.tempdir,"structure.lmpstruct"))
    shutil.copyfile(os.path.join(_getLAMMPSResourceDirectory(), "calc_energy.lmpin"), os.path.join(self.tempdir,"calc_energy.lmpin"))

    inputExpect = [
      (1.24246478E-02, "Al Al"), # Al Al
      (-0.121863537, "Al Cu"),  # Al Cu
      (-0.179150283, "Cu Cu") # Cu Cu
    ]

    oldpwd = os.getcwd()
    os.chdir(self.tempdir)
    try:
      eampotentials, pairPotentials = exampleModule.makePotentialObjects()
      eamPotentials = None

      # Create EAMPotential
      def density(r):
        return 0.0

      def embed(rho):
        return 0.0

      eamPotentials = [
        EAMPotential("Cu", 29, 63.55, embed, density),
        EAMPotential("Al", 13, 26.98, embed, density)]

      nrho = 50000
      drho = 0.001

      nr = 12000
      dr = 0.001

      from atsim.potentials import writeSetFL
      with open("table.set", 'wb') as outfile:
        writeSetFL(
          nrho, drho,
          nr, dr,
          eamPotentials,
          pairPotentials,
          comments = ['Zhou Al Cu', "", ""],
          out= outfile)

      for expect, potmap in inputExpect:
        with open("potentials.lmpinc", "wb") as potfile:
          print >>potfile, "pair_style eam/alloy"
          print >>potfile, "pair_coeff * * table.set "+potmap
        runLAMMPS()
        energy = extractLAMMPSEnergy()
        self.assertAlmostEquals(expect, energy, msg = potmap)
    finally:
      os.chdir(oldpwd)


  @needsLAMMPS
  def testExampleA_Density(self):
    """Test density functions correctly defined for EAM tabulation documentation example 2a"""
    exampleModule = _loadModule(self.test_nameA)

    shutil.copyfile(os.path.join(_getLAMMPSResourceDirectory(), "setfl_pair.lmpstruct"), os.path.join(self.tempdir,"structure.lmpstruct"))
    shutil.copyfile(os.path.join(_getLAMMPSResourceDirectory(), "calc_energy.lmpin"), os.path.join(self.tempdir,"calc_energy.lmpin"))

    oldpwd = os.getcwd()
    os.chdir(self.tempdir)
    try:
      def nullfunc(r):
        return 0.0

      eamPotentials, pairPotentials = exampleModule.makePotentialObjects()
      pairPotentials = None

      pairPotentials = [
        Potential("Cu", "Cu", nullfunc),
        Potential("Al", "Al", nullfunc),
        Potential("Al", "Cu", nullfunc)
      ]

      # Create EAMPotential
      def embed(rho):
        return rho

      for epot in eamPotentials:
        epot.embeddingFunction = embed

      nrho = 50000
      drho = 0.001

      nr = 12000
      dr = 0.001

      from atsim.potentials import writeSetFL

      with open("table.set", 'wb') as outfile:
        writeSetFL(
          nrho, drho,
          nr, dr,
          eamPotentials,
          pairPotentials,
          comments = ['Zhou Al Cu', "", ""],
          out= outfile)
      # import pdb;pdb.set_trace()

      inputExpect = [
        ( 2.0*1.218017211, "Cu Cu"),  # Cu Cu
        ( 2.0*1.716990097, "Al Al"), # Al Al
        ( 1.218017211+1.716990097, "Al Cu")  # Al Cu
      ]

      for expect, potmap in inputExpect:
        with open("potentials.lmpinc", "wb") as potfile:
          print >>potfile, "pair_style eam/alloy"
          print >>potfile, "pair_coeff * * table.set "+potmap
        runLAMMPS()
        energy = extractLAMMPSEnergy()
        self.assertAlmostEquals(expect, energy, msg = potmap)

      # Now repeat for triplet of atoms
      shutil.copyfile(os.path.join(_getLAMMPSResourceDirectory(), "setfl_triplet.lmpstruct"), os.path.join(self.tempdir,"structure.lmpstruct"))
      dens_Cu = [p.electronDensityFunction for p in eamPotentials if p.species == "Cu"][0]
      dens_Al = [p.electronDensityFunction for p in eamPotentials if p.species == "Al"][0]
      hyp = 3.818376618407357
      inputExpect = [
        ( 4* dens_Cu(2.7) + 2* dens_Cu(hyp), "Cu Cu"),  # Cu Cu
        ( 4* dens_Al(2.7) + 2* dens_Al(hyp), "Al Al"), # Al Al
        ( 2* dens_Cu(2.7) + 2* dens_Cu(hyp) + 2*dens_Al(2.7), "Al Cu"),  # Al Cu
        ( 2* dens_Al(2.7) + 2* dens_Al(hyp) + 2*dens_Cu(2.7), "Cu Al")  # Cu Al
      ]
      for expect, potmap in inputExpect:
        with open("potentials.lmpinc", "wb") as potfile:
          print >>potfile, "pair_style eam/alloy"
          print >>potfile, "pair_coeff * * table.set "+potmap
        runLAMMPS()
        energy = extractLAMMPSEnergy()
        self.assertAlmostEquals(expect, energy, msg = potmap)

    finally:
      os.chdir(oldpwd)


  def testWriteSetFLEmbedCu(self):
    """Test Cu embedding function correctly defined for EAM tabulation documentation example 2"""
    exampleModule = _loadModule(self.test_nameA)

    eamPotentials, pairPotentials = exampleModule.makePotentialObjects()
    embed_Cu = [p.embeddingFunction for p in eamPotentials if p.species == "Cu"][0]

    self.assertAlmostEquals(-1.76619128240398, embed_Cu(10.0))
    self.assertAlmostEquals(-2.18790796129658, embed_Cu(20.0))
    self.assertAlmostEquals(-2.17281697911785, embed_Cu(30.0))
    self.assertAlmostEquals(-2.13787794765212, embed_Cu(40.0))


  def testWriteSetFLEmbedAl(self):
    """Test Al embedding function correctly defined for EAM tabulation documentation example 2"""
    exampleModule = _loadModule(self.test_nameA)

    eamPotentials, pairPotentials = exampleModule.makePotentialObjects()
    embed_Al = [p.embeddingFunction for p in eamPotentials if p.species == "Al"][0]

    self.assertAlmostEquals(-2.35881750559297, embed_Al(10.0))
    self.assertAlmostEquals(-2.82971758138417, embed_Al(20.0))
    self.assertAlmostEquals(-2.75841139984064, embed_Al(30.0))
    self.assertAlmostEquals(-2.47821972143384, embed_Al(40.0))


  @needsLAMMPS
  def testExampleA(self):
    """Test example eam_tabulate_example2a.py"""
    exampleModule = _loadModule(self.test_nameA)

    oldpwd = os.getcwd()
    os.chdir(self.tempdir)
    try:
      shutil.copyfile(os.path.join(_getLAMMPSResourceDirectory(), "calc_energy.lmpin"), os.path.join(self.tempdir,"calc_energy.lmpin"))
      shutil.copyfile(os.path.join(_getLAMMPSResourceDirectory(), "setfl_triplet.lmpstruct"), os.path.join(self.tempdir,"structure.lmpstruct"))
      shutil.copyfile(os.path.join(_getLAMMPSResourceDirectory(), "Zhou_AlCu.setfl"), os.path.join(self.tempdir, "table.setfl"))

      potmaps = ["Cu Cu","Al Al", "Al Cu", "Cu Al"]
      expect = []

      # Run the Zhou tabulation created using tools from http://www.ctcms.nist.gov/potentials/Zhou04.html
      for potmap in potmaps:
        with open("potentials.lmpinc", "wb") as potfile:
          print >>potfile, "pair_style eam/alloy"
          print >>potfile, "pair_coeff * * table.setfl "+potmap
        runLAMMPS()
        energy = extractLAMMPSEnergy()
        self.assertTrue(energy != None)
        expect.append(energy)

      # Make  tabulation
      nrho = 50000
      drho = 0.001

      nr = 12000
      dr = 0.001

      # from atsim.potentials import writeSetFL
      # pairPotentials = self.createSetFlPairPots()
      # eamPotentials  = self.createSetFLEAMPots()

      # with open("table.set", 'wb') as outfile:
      #   writeSetFL(
      #     nrho, drho,
      #     nr, dr,
      #     eamPotentials,
      #     pairPotentials,
      #     comments = ['Zhou Al Cu', "", ""],
      #     out= outfile)

      exampleModule.main()

      for potmap,expectEnergy in zip(potmaps,expect):
        with open("potentials.lmpinc", "wb") as potfile:
          print >>potfile, "pair_style eam/alloy"
          print >>potfile, "pair_coeff * * Zhou_AlCu.setfl "+potmap
        runLAMMPS()
        energy = extractLAMMPSEnergy()
        self.assertAlmostEquals(expectEnergy, energy, places = 4, msg = potmap)

    finally:
      os.chdir(oldpwd)

  @needsLAMMPS
  @needsDLPOLY
  def testCrossCheckLAMMPS_DLPOLY(self):
    """Check that models tabulated for LAMMPS and DL_POLY give the same result (cross check example 2a and 2b)."""
    exampleModuleA = _loadModule(self.test_nameA)
    exampleModuleB = _loadModule(self.test_nameB)

    oldpwd = os.getcwd()
    os.chdir(self.tempdir)
    try:
      # DL_POLY Tabulation
      shutil.copyfile(os.path.join(_getDLPolyResourceDirectory(), "CONTROL_random_Al_Cu"), os.path.join(self.tempdir,"CONTROL"))
      shutil.copyfile(os.path.join(_getDLPolyResourceDirectory(), "CONFIG_random_Al_Cu"), os.path.join(self.tempdir,"CONFIG"))
      shutil.copyfile(os.path.join(_getDLPolyResourceDirectory(), "FIELD_random_Al_Cu"), os.path.join(self.tempdir,"FIELD"))

      # Create TABEAM
      exampleModuleB.main()

      runDLPoly()
      # import pdb;pdb.set_trace()
      dlpolyenergy = extractDLPOLYEnergy()

      # LAMMPS Tabulation
      shutil.copyfile(os.path.join(_getLAMMPSResourceDirectory(), "calc_energy.lmpin"), os.path.join(self.tempdir,"calc_energy.lmpin"))
      shutil.copyfile(os.path.join(_getLAMMPSResourceDirectory(), "random_Al_Cu.lmpstruct"), os.path.join(self.tempdir,"structure.lmpstruct"))
      with open("potentials.lmpinc", "wb") as potfile:
        print >>potfile, "pair_style eam/alloy"
        print >>potfile, "pair_coeff * * Zhou_AlCu.setfl Al Cu"

      # Create the table files
      exampleModuleA.main()

      # import pdb;pdb.set_trace()
      runLAMMPS()
      lammpsenergy = extractLAMMPSEnergy()

      self.assertAlmostEquals(lammpsenergy, dlpolyenergy, places = 4)

    finally:
      os.chdir(oldpwd)



#eam_tabulate_example3a.py
class eam_tabulate_example3TestCase(TempfileTestCase):
  """Test docs/potentials/eam_tabulate_example3a.py and eam_tabulate_example3b.py"""
  test_nameA = os.path.join(_getDocsDirectory(), "eam_tabulate_example3a.py")
  test_nameB = os.path.join(_getDocsDirectory(), "eam_tabulate_example3b.py")

  def testExampleA(self):
    """Test example eam_tabulate_example3a.py"""
    exampleModule = _loadModule(self.test_nameA)

    oldpwd = os.getcwd()
    os.chdir(self.tempdir)
    try:
      exampleModule.main()
    finally:
      os.chdir(oldpwd)


  def testExampleB(self):
    """Test example eam_tabulate_example3b.py"""
    exampleModule = _loadModule(self.test_nameB)

    oldpwd = os.getcwd()
    os.chdir(self.tempdir)
    try:
      exampleModule.main()
    finally:
      os.chdir(oldpwd)

  @needsLAMMPS
  @needsDLPOLY
  def testCrossCheckLAMMPS_DLPOLY(self):
    """Check that models tabulated for LAMMPS and DL_POLY give the same result"""
    exampleModuleA = _loadModule(self.test_nameA)
    exampleModuleB = _loadModule(self.test_nameB)

    oldpwd = os.getcwd()
    os.chdir(self.tempdir)
    try:
      # DL_POLY Tabulation
      shutil.copyfile(os.path.join(_getDLPolyResourceDirectory(), "CONTROL_random_Al_Fe"), os.path.join(self.tempdir,"CONTROL"))
      shutil.copyfile(os.path.join(_getDLPolyResourceDirectory(), "CONFIG_random_Al_Fe"), os.path.join(self.tempdir,"CONFIG"))
      shutil.copyfile(os.path.join(_getDLPolyResourceDirectory(), "FIELD_random_Al_Fe"), os.path.join(self.tempdir,"FIELD"))

      # Create TABEAM
      exampleModuleB.main()

      # import pdb;pdb.set_trace()
      runDLPoly()
      dlpolyenergy = extractDLPOLYEnergy()

      # LAMMPS Tabulation
      shutil.copyfile(os.path.join(_getLAMMPSResourceDirectory(), "calc_energy.lmpin"), os.path.join(self.tempdir,"calc_energy.lmpin"))
      shutil.copyfile(os.path.join(_getLAMMPSResourceDirectory(), "random_Al_Fe.lmpstruct"), os.path.join(self.tempdir,"structure.lmpstruct"))
      with open("potentials.lmpinc", "wb") as potfile:
        print >>potfile, "pair_style eam/fs"
        print >>potfile, "pair_coeff * * Mendelev_Al_Fe.eam.fs Al Fe"

      # Create the table files
      exampleModuleA.main()

      # import pdb;pdb.set_trace()
      runLAMMPS()
      lammpsenergy = extractLAMMPSEnergy()

      self.assertAlmostEquals(lammpsenergy, dlpolyenergy, places = 4)

    finally:
      os.chdir(oldpwd)

try:
  import numpy
  NUMPY_AVAILABLE = True
except ImportError:
  NUMPY_AVAILABLE = False


#zbl_spline.py
class zbl_splineTestCase(TempfileTestCase):
  """Test docs/potentials/zbl_spline.py"""
  test_name = os.path.join(_getDocsDirectory(), "zbl_spline.py")

  @unittest.skipIf(not NUMPY_AVAILABLE, "numpy is not installed")
  def testExample(self):
    """Test example zbl_spline.py"""
    exampleModule = _loadModule(self.test_name)
    oldpwd = os.getcwd()
    os.chdir(self.tempdir)
    try:
      exampleModule.main()
    finally:
      os.chdir(oldpwd)

