"""Tests to make sure tests given in the documentation work"""
from __future__ import absolute_import, print_function

import os
import shutil
import sys
import unittest
from builtins import range, str, zip

import atsim.potentials
import atsim.potentials.config
from atsim.potentials import EAMPotential, Potential

import py.path
import pytest

from ._rundlpoly import extractDLPOLYEnergy, needsDLPOLY, runDLPoly
from ._rungulp import gulp_uo2_energy_fixture, needsGULP
from ._runlammps import extractLAMMPSEnergy, needsLAMMPS, runLAMMPS
from ._tempfiletestcase import TempfileTestCase


def _getDocsDirectory():
    """Returns absolute path to docs/ directory"""
    docsdir = os.path.join("docs", "potentials")
    return os.path.abspath(docsdir)


def _get_user_guide_directory():
    """Returns absolute path to docs/ directory"""
    docsdir = os.path.join("docs", "user_guide")
    return os.path.abspath(docsdir)


def _getLAMMPSResourceDirectory():
    return os.path.join(os.path.dirname(__file__), 'lammps_resources')


def _getDLPolyResourceDirectory():
    return os.path.join(os.path.dirname(__file__), 'dl_poly_resources')


if sys.version_info.major == 3 and sys.version_info.minor >= 5:

    import importlib.util

    def _loadModule(scriptname):
        name = os.path.basename(scriptname)
        name = os.path.splitext(name)[0]
        spec = importlib.util.spec_from_file_location(name, scriptname)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod

elif sys.version_info.major == 2:

    import imp

    def _loadModule(scriptname):
        name = os.path.basename(scriptname)
        name = os.path.splitext(name)[0]
        with open(scriptname) as infile:
            mod = imp.load_module(name, infile, scriptname, ('.py', 'U', 1))
            # mod = importlib.load_module(name,infile, scriptname, ('.py', 'U', 1))
        return mod

else:
    raise Exception(
        "No implementation of _loadModule for python version {}".format(sys.version))

# basak_tabulate.py


class basak_tabulateTestCase(TempfileTestCase):
    """Test docs/potentials/basak_tabulate.py"""
    test_name = os.path.join(
        _getDocsDirectory(), os.pardir, "user_guide", "basak_tabulate.py")

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
                spairs = [
                    tuple(sorted(p)) for p in itertools.combinations_with_replacement(species, 2)]
                species = set(spairs)

                pobjs = [p for p in pobjs if tuple(
                    sorted([p.speciesA, p.speciesB])) in species]

                with open('TABLE', 'w') as outfile:
                    atsim.potentials.writePotentials(
                        'DL_POLY',
                        pobjs,
                        6.5, 6500,
                        out=outfile)

                d = dict(d)
                for i in range(4):
                    r = float(i)
                    r += 1.0
                    expect = testfunc(r)

                    d['Ax'] = r

                    with open(os.path.join(dldir, "CONFIG_pair.in"), 'r') as infile:
                        with open("CONFIG", "w") as outfile:
                            outfile.write(infile.read() % d)

                    from io import StringIO

                    sio = StringIO()
                    print(u"vdw %d" % len(pobjs), file=sio)
                    for p in pobjs:
                        print(u"%s %s tab" %
                              (p.speciesA, p.speciesB), file=sio)
                    d["potDef"] = sio.getvalue()

                    with open(os.path.join(dldir, "FIELD_pair.in"), 'r') as infile:
                        with open("FIELD", "w") as outfile:
                            outfile.write(infile.read() % d)

                    runDLPoly()
                    dlenergy = extractDLPOLYEnergy()
                    self.assertAlmostEqual(expect, dlenergy, places=4)
                    os.remove("STATIS")
                    os.remove("CONFIG")
                    os.remove("FIELD")

            d = dict(speciesA="O", speciesB="O",
                     Ax=0.0, Ay=0.0, Az=0.0,
                     Bx=0.0, By=0.0, Bz=0.0)
            testfunc = potentialforms.buck(1633.00510, 0.327022, 3.948790)
            dotest(d, testfunc)

            d = dict(speciesA="U", speciesB="U",
                     Ax=0.0, Ay=0.0, Az=0.0,
                     Bx=0.0, By=0.0, Bz=0.0)
            testfunc = potentialforms.buck(294.640000, 0.327022, 0.0)
            dotest(d, testfunc)

            d = dict(speciesA="O", speciesB="U",
                     Ax=0.0, Ay=0.0, Az=0.0,
                     Bx=0.0, By=0.0, Bz=0.0)
            buck_OU = potentialforms.buck(693.648700, 0.327022, 0.0)
            morse_OU = potentialforms.morse(1.6500, 2.36900, 0.577190)
            testfunc = atsim.potentials.plus(buck_OU, morse_OU)

            dotest(d, testfunc)

        finally:
            os.chdir(oldpwd)

    def testLAMMPSExample(self):
        """Test doumentation example Quick-Start: LAMMPS."""

        test_name = os.path.join(
            _get_user_guide_directory(), "basak_tabulate_lammps.py")

        oldpwd = os.getcwd()
        try:
            os.chdir(self.tempdir)
            exampleModule = _loadModule(test_name)
            exampleModule.main()
        finally:
            os.chdir(oldpwd)

# eam_tabulate_example1.py


class eam_tabulate_example1TestCase(TempfileTestCase):
    """Test docs/potentials/eam_example1.py"""
    test_name = os.path.join(_getDocsDirectory(), "eam_example1.py")

    @needsLAMMPS
    def testExample(self):
        """Test example eam_example1.py"""
        exampleModule = _loadModule(self.test_name)

        shutil.copyfile(os.path.join(_getLAMMPSResourceDirectory(
        ), "doc_example1_Ag_fcc.lmpstruct"), os.path.join(self.tempdir, "structure.lmpstruct"))
        shutil.copyfile(os.path.join(_getLAMMPSResourceDirectory(
        ), "calc_energy.lmpin"), os.path.join(self.tempdir, "calc_energy.lmpin"))

        oldpwd = os.getcwd()
        os.chdir(self.tempdir)
        try:
            with open("potentials.lmpinc", "w") as potfile:
                potfile.write("pair_style eam/alloy\n")
                potfile.write("pair_coeff * * Ag.eam.alloy Ag\n")

            # Run the main method
            exampleModule.main()

            runLAMMPS()
            energy = extractLAMMPSEnergy()

            self.assertAlmostEqual(-8.23982879, energy, places=5)
        finally:
            os.chdir(oldpwd)


class eam_tabulate_example1_procedural_TestCase(TempfileTestCase):
    """Test docs/potentials/eam_tabulate_example1.py"""
    test_name = os.path.join(_getDocsDirectory(), "eam_tabulate_example1.py")

    @needsLAMMPS
    def testExample(self):
        """Test example eam_tabulate_example1.py"""
        exampleModule = _loadModule(self.test_name)

        shutil.copyfile(os.path.join(_getLAMMPSResourceDirectory(
        ), "doc_example1_Ag_fcc.lmpstruct"), os.path.join(self.tempdir, "structure.lmpstruct"))
        shutil.copyfile(os.path.join(_getLAMMPSResourceDirectory(
        ), "calc_energy.lmpin"), os.path.join(self.tempdir, "calc_energy.lmpin"))

        oldpwd = os.getcwd()
        os.chdir(self.tempdir)
        try:
            with open("potentials.lmpinc", "w") as potfile:
                potfile.write("pair_style eam\n")
                potfile.write("pair_coeff 1 1 Ag.eam\n")

            # Run the main method
            exampleModule.main()

            runLAMMPS()
            energy = extractLAMMPSEnergy()

            self.assertAlmostEqual(-8.23982879, energy, places=5)
        finally:
            os.chdir(oldpwd)


class eam_tabulate_example2TestCase(TempfileTestCase):
    """Test docs/potentials/eam_tabulate_example2a.py and docs/potentials/eam_tabulate_example2b.py"""
    test_nameA = os.path.join(_getDocsDirectory(), "eam_tabulate_example2a.py")
    test_nameA_obj = os.path.join(_getDocsDirectory(), "eam_example2a.py")
    test_nameB = os.path.join(_getDocsDirectory(), "eam_tabulate_example2b.py")
    test_nameB_obj = os.path.join(_getDocsDirectory(), "eam_example2b.py")

    @needsLAMMPS
    def testExampleA_Pair(self):
        """Test pair-potentials correctly defined for EAM tabulation documentation example 2a"""
        exampleModule = _loadModule(self.test_nameA)

        shutil.copyfile(os.path.join(_getLAMMPSResourceDirectory(
        ), "setfl_pair.lmpstruct"), os.path.join(self.tempdir, "structure.lmpstruct"))
        shutil.copyfile(os.path.join(_getLAMMPSResourceDirectory(
        ), "calc_energy.lmpin"), os.path.join(self.tempdir, "calc_energy.lmpin"))

        inputExpect = [
            (1.24246478E-02, "Al Al"),  # Al Al
            (-0.121863537, "Al Cu"),  # Al Cu
            (-0.179150283, "Cu Cu")  # Cu Cu
        ]

        oldpwd = os.getcwd()
        os.chdir(self.tempdir)
        try:
            _eampotentials, pairPotentials = exampleModule.makePotentialObjects()
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
            with open("table.set", 'w') as outfile:
                writeSetFL(
                    nrho, drho,
                    nr, dr,
                    eamPotentials,
                    pairPotentials,
                    comments=['Zhou Al Cu', "", ""],
                    out=outfile)

            for expect, potmap in inputExpect:
                with open("potentials.lmpinc", "w") as potfile:
                    potfile.write("pair_style eam/alloy\n")
                    potfile.write("pair_coeff * * table.set "+potmap+"\n")
                runLAMMPS()
                energy = extractLAMMPSEnergy()
                self.assertAlmostEqual(expect, energy, msg=potmap)
        finally:
            os.chdir(oldpwd)

    @needsLAMMPS
    def testExampleA_Density(self):
        """Test density functions correctly defined for EAM tabulation documentation example 2a"""
        exampleModule = _loadModule(self.test_nameA)

        shutil.copyfile(os.path.join(_getLAMMPSResourceDirectory(
        ), "setfl_pair.lmpstruct"), os.path.join(self.tempdir, "structure.lmpstruct"))
        shutil.copyfile(os.path.join(_getLAMMPSResourceDirectory(
        ), "calc_energy.lmpin"), os.path.join(self.tempdir, "calc_energy.lmpin"))

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

            with open("table.set", 'w') as outfile:
                writeSetFL(
                    nrho, drho,
                    nr, dr,
                    eamPotentials,
                    pairPotentials,
                    comments=['Zhou Al Cu', "", ""],
                    out=outfile)

            inputExpect = [
                (2.0*1.218017211, "Cu Cu"),  # Cu Cu
                (2.0*1.716990097, "Al Al"),  # Al Al
                (1.218017211+1.716990097, "Al Cu")  # Al Cu
            ]

            for expect, potmap in inputExpect:
                with open("potentials.lmpinc", "w") as potfile:
                    potfile.write("pair_style eam/alloy\n")
                    potfile.write("pair_coeff * * table.set "+potmap+"\n")
                runLAMMPS()
                energy = extractLAMMPSEnergy()
                self.assertAlmostEqual(expect, energy, msg=potmap)

            # Now repeat for triplet of atoms
            shutil.copyfile(os.path.join(_getLAMMPSResourceDirectory(
            ), "setfl_triplet.lmpstruct"), os.path.join(self.tempdir, "structure.lmpstruct"))
            dens_Cu = [
                p.electronDensityFunction for p in eamPotentials if p.species == "Cu"][0]
            dens_Al = [
                p.electronDensityFunction for p in eamPotentials if p.species == "Al"][0]
            hyp = 3.818376618407357
            inputExpect = [
                (4 * dens_Cu(2.7) + 2 * dens_Cu(hyp), "Cu Cu"),  # Cu Cu
                (4 * dens_Al(2.7) + 2 * dens_Al(hyp), "Al Al"),  # Al Al
                (2 * dens_Cu(2.7) + 2 * dens_Cu(hyp) + \
                    2*dens_Al(2.7), "Al Cu"),  # Al Cu
                (2 * dens_Al(2.7) + 2 * dens_Al(hyp) + \
                    2*dens_Cu(2.7), "Cu Al")  # Cu Al
            ]
            for expect, potmap in inputExpect:
                with open("potentials.lmpinc", "w") as potfile:
                    potfile.write("pair_style eam/alloy\n")
                    potfile.write("pair_coeff * * table.set "+potmap+"\n")
                runLAMMPS()
                energy = extractLAMMPSEnergy()
                self.assertAlmostEqual(expect, energy, msg=potmap)

        finally:
            os.chdir(oldpwd)

    def testWriteSetFLEmbedCu(self):
        """Test Cu embedding function correctly defined for EAM tabulation documentation example 2"""
        exampleModule = _loadModule(self.test_nameA)

        eamPotentials, _pairPotentials = exampleModule.makePotentialObjects()
        embed_Cu = [
            p.embeddingFunction for p in eamPotentials if p.species == "Cu"][0]

        self.assertAlmostEqual(-1.76619128240398, embed_Cu(10.0))
        self.assertAlmostEqual(-2.18790796129658, embed_Cu(20.0))
        self.assertAlmostEqual(-2.17281697911785, embed_Cu(30.0))
        self.assertAlmostEqual(-2.13787794765212, embed_Cu(40.0))

    def testWriteSetFLEmbedAl(self):
        """Test Al embedding function correctly defined for EAM tabulation documentation example 2"""
        exampleModule = _loadModule(self.test_nameA)

        eamPotentials, _pairPotentials = exampleModule.makePotentialObjects()
        embed_Al = [
            p.embeddingFunction for p in eamPotentials if p.species == "Al"][0]

        self.assertAlmostEqual(-2.35881750559297, embed_Al(10.0))
        self.assertAlmostEqual(-2.82971758138417, embed_Al(20.0))
        self.assertAlmostEqual(-2.75841139984064, embed_Al(30.0))
        self.assertAlmostEqual(-2.47821972143384, embed_Al(40.0))

    @needsLAMMPS
    def testExampleA(self):
        """Test example eam_tabulate_example2a.py"""
        exampleModule = _loadModule(self.test_nameA)

        oldpwd = os.getcwd()
        os.chdir(self.tempdir)
        try:
            shutil.copyfile(os.path.join(_getLAMMPSResourceDirectory(
            ), "calc_energy.lmpin"), os.path.join(self.tempdir, "calc_energy.lmpin"))
            shutil.copyfile(os.path.join(_getLAMMPSResourceDirectory(
            ), "setfl_triplet.lmpstruct"), os.path.join(self.tempdir, "structure.lmpstruct"))
            shutil.copyfile(os.path.join(_getLAMMPSResourceDirectory(
            ), "Zhou_AlCu.setfl"), os.path.join(self.tempdir, "table.setfl"))

            potmaps = ["Cu Cu", "Al Al", "Al Cu", "Cu Al"]
            expect = []

            # Run the Zhou tabulation created using tools from http://www.ctcms.nist.gov/potentials/Zhou04.html
            for potmap in potmaps:
                with open("potentials.lmpinc", "w") as potfile:
                    potfile.write("pair_style eam/alloy\n")
                    potfile.write("pair_coeff * * table.setfl "+potmap+"\n")
                runLAMMPS()
                energy = extractLAMMPSEnergy()
                self.assertTrue(energy != None)
                expect.append(energy)

            # Make  tabulation
            exampleModule.main()

            for potmap, expectEnergy in zip(potmaps, expect):
                with open("potentials.lmpinc", "w") as potfile:
                    potfile.write("pair_style eam/alloy\n")
                    potfile.write(
                        "pair_coeff * * Zhou_AlCu.eam.alloy "+potmap+"\n")
                runLAMMPS()
                energy = extractLAMMPSEnergy()
                self.assertAlmostEqual(
                    expectEnergy, energy, places=4, msg=potmap)

        finally:
            os.chdir(oldpwd)

    def testExampleA_obj(self):
        """Test example eam_example2a.py"""
        exampleModule = _loadModule(self.test_nameA_obj)

        oldpwd = os.getcwd()
        os.chdir(self.tempdir)
        try:
            shutil.copyfile(os.path.join(_getLAMMPSResourceDirectory(
            ), "calc_energy.lmpin"), os.path.join(self.tempdir, "calc_energy.lmpin"))
            shutil.copyfile(os.path.join(_getLAMMPSResourceDirectory(
            ), "setfl_triplet.lmpstruct"), os.path.join(self.tempdir, "structure.lmpstruct"))
            shutil.copyfile(os.path.join(_getLAMMPSResourceDirectory(
            ), "Zhou_AlCu.setfl"), os.path.join(self.tempdir, "table.setfl"))

            potmaps = ["Cu Cu", "Al Al", "Al Cu", "Cu Al"]
            expect = []

            # Run the Zhou tabulation created using tools from http://www.ctcms.nist.gov/potentials/Zhou04.html
            for potmap in potmaps:
                with open("potentials.lmpinc", "w") as potfile:
                    potfile.write("pair_style eam/alloy\n")
                    potfile.write("pair_coeff * * table.setfl "+potmap+"\n")
                runLAMMPS()
                energy = extractLAMMPSEnergy()
                self.assertTrue(energy != None)
                expect.append(energy)

            exampleModule.main()

            for potmap, expectEnergy in zip(potmaps, expect):
                with open("potentials.lmpinc", "w") as potfile:
                    potfile.write("pair_style eam/alloy\n")
                    potfile.write(
                        "pair_coeff * * Zhou_AlCu.eam.alloy "+potmap+"\n")
                runLAMMPS()
                energy = extractLAMMPSEnergy()
                self.assertAlmostEqual(
                    expectEnergy, energy, places=4, msg=potmap)

        finally:
            os.chdir(oldpwd)

    @needsLAMMPS
    @needsDLPOLY
    def testExample2b(self):
        exampleModuleB = _loadModule(self.test_nameB)
        self._crossCheckLAMMPS_DLPOLY(exampleModuleB)

    @needsLAMMPS
    @needsDLPOLY
    def testExample2b_obj(self):
        exampleModuleB = _loadModule(self.test_nameB_obj)
        self._crossCheckLAMMPS_DLPOLY(exampleModuleB)

    def _crossCheckLAMMPS_DLPOLY(self, exampleModuleB):
        """Check that models tabulated for LAMMPS and DL_POLY give the same result (cross check example 2a and 2b)."""
        exampleModuleA = _loadModule(self.test_nameA)
        # exampleModuleB = _loadModule(self.test_nameB)

        oldpwd = os.getcwd()
        os.chdir(self.tempdir)
        try:
            # DL_POLY Tabulation
            shutil.copyfile(os.path.join(_getDLPolyResourceDirectory(
            ), "CONTROL_random_Al_Cu"), os.path.join(self.tempdir, "CONTROL"))
            shutil.copyfile(os.path.join(_getDLPolyResourceDirectory(
            ), "CONFIG_random_Al_Cu"), os.path.join(self.tempdir, "CONFIG"))
            shutil.copyfile(os.path.join(_getDLPolyResourceDirectory(
            ), "FIELD_random_Al_Cu"), os.path.join(self.tempdir, "FIELD"))

            # Create TABEAM
            exampleModuleB.main()

            runDLPoly()
            # import pdb;pdb.set_trace()
            dlpolyenergy = extractDLPOLYEnergy()

            # LAMMPS Tabulation
            shutil.copyfile(os.path.join(_getLAMMPSResourceDirectory(
            ), "calc_energy.lmpin"), os.path.join(self.tempdir, "calc_energy.lmpin"))
            shutil.copyfile(os.path.join(_getLAMMPSResourceDirectory(
            ), "random_Al_Cu.lmpstruct"), os.path.join(self.tempdir, "structure.lmpstruct"))
            with open("potentials.lmpinc", "w") as potfile:
                potfile.write("pair_style eam/alloy\n")
                potfile.write("pair_coeff * * Zhou_AlCu.eam.alloy Al Cu\n")

            # Create the table files
            exampleModuleA.main()

            # import pdb;pdb.set_trace()
            runLAMMPS()
            lammpsenergy = extractLAMMPSEnergy()

            self.assertAlmostEqual(lammpsenergy, dlpolyenergy, places=4)

        finally:
            os.chdir(oldpwd)

# eam_tabulate_example3a.py
# class eam_tabulate_example3TestCase(TempfileTestCase):
#   """Test docs/potentials/eam_tabulate_example3a.py and eam_tabulate_example3b.py"""
#   test_nameA = os.path.join(_getDocsDirectory(), "eam_tabulate_example3a.py")
#   test_nameA_obj = os.path.join(_getDocsDirectory(), "eam_example3a.py")
#   test_nameB = os.path.join(_getDocsDirectory(), "eam_tabulate_example3b.py")
#   test_nameB_obj = os.path.join(_getDocsDirectory(), "eam_example3b.py")


def test_eam_tabulate_example3_ExampleA(tmpdir):
    """Test example eam_tabulate_example3a.py"""
    test_nameA = os.path.join(_getDocsDirectory(), "eam_tabulate_example3a.py")
    exampleModule = _loadModule(test_nameA)

    oldpwd = os.getcwd()
    os.chdir(tmpdir.strpath)
    try:
        exampleModule.main()
    finally:
        os.chdir(oldpwd)


def test_eam_tabulate_example3_ExampleB(tmpdir):
    """Test example eam_tabulate_example3b.py"""
    test_nameB = os.path.join(_getDocsDirectory(), "eam_tabulate_example3b.py")
    exampleModule = _loadModule(test_nameB)

    oldpwd = os.getcwd()
    os.chdir(tmpdir.strpath)
    try:
        exampleModule.main()
    finally:
        os.chdir(oldpwd)


@needsLAMMPS
@needsDLPOLY
@pytest.mark.parametrize(("test_name_A", "test_name_B"), [
    ("eam_tabulate_example3a.py", "eam_tabulate_example3b.py"),
    ("eam_example3a.py", "eam_example3b.py"),
    ("eam_example3a.py", "eam_tabulate_example3b.py"),
    ("eam_tabulate_example3a.py", "eam_example3b.py")
])
def test_eam_tabulate_example3_cross_check_LAMMPS_DLPOLY(test_name_A, test_name_B, tmpdir):
    test_name_A = os.path.join(_getDocsDirectory(), test_name_A)
    test_name_B = os.path.join(_getDocsDirectory(), test_name_B)

    exampleModuleA = _loadModule(test_name_A)
    exampleModuleB = _loadModule(test_name_B)
    _crossCheckLAMMPS_DLPOLY(tmpdir, exampleModuleA, exampleModuleB)


def _crossCheckLAMMPS_DLPOLY(tmpdir, exampleModuleA, exampleModuleB):
    """Check that models tabulated for LAMMPS and DL_POLY give the same result"""
    tmpdir = tmpdir.strpath
    oldpwd = os.getcwd()
    os.chdir(tmpdir)
    try:
        # DL_POLY Tabulation
        shutil.copyfile(os.path.join(_getDLPolyResourceDirectory(
        ), "CONTROL_random_Al_Fe"), os.path.join(tmpdir, "CONTROL"))
        shutil.copyfile(os.path.join(_getDLPolyResourceDirectory(
        ), "CONFIG_random_Al_Fe"), os.path.join(tmpdir, "CONFIG"))
        shutil.copyfile(os.path.join(_getDLPolyResourceDirectory(
        ), "FIELD_random_Al_Fe"), os.path.join(tmpdir, "FIELD"))

        # Create TABEAM
        exampleModuleB.main()

        runDLPoly()
        dlpolyenergy = extractDLPOLYEnergy()

        # LAMMPS Tabulation
        shutil.copyfile(os.path.join(_getLAMMPSResourceDirectory(
        ), "calc_energy.lmpin"), os.path.join(tmpdir, "calc_energy.lmpin"))
        shutil.copyfile(os.path.join(_getLAMMPSResourceDirectory(
        ), "random_Al_Fe.lmpstruct"), os.path.join(tmpdir, "structure.lmpstruct"))
        with open("potentials.lmpinc", "w") as potfile:
            potfile.write("pair_style eam/fs\n")
            potfile.write("pair_coeff * * Mendelev_Al_Fe.eam.fs Al Fe\n")

        # Create the table files
        exampleModuleA.main()

        # import pdb;pdb.set_trace()
        runLAMMPS()
        lammpsenergy = extractLAMMPSEnergy()

        assert lammpsenergy == pytest.approx(dlpolyenergy)

    finally:
        os.chdir(oldpwd)


try:
    import numpy
    NUMPY_AVAILABLE = True
except ImportError:
    NUMPY_AVAILABLE = False

# zbl_spline.py


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

            output_path = os.path.join(self.tempdir, "bks_buck.dat")
            assert os.path.exists(output_path)

            with open(output_path) as infile:
                line = next(infile)
                tokens = line.split()
                r = float(tokens[0])
                assert pytest.approx(0.1) == r

                line = next(infile)
                tokens = line.split()
                r = float(tokens[0])
                assert pytest.approx(0.1 + (10.0-0.1)/5000.0) == r

        finally:
            os.chdir(oldpwd)


basak_aspot_files = py.path.local(__file__).dirpath(
    "..", "docs", "user_guide", "example_files").listdir("basak*.aspot")


@needsLAMMPS
@pytest.mark.parametrize("aspotfile", basak_aspot_files)
def test_basak_files(tmpdir, aspotfile):
    # Copy files in from example directory
    srcdir = py.path.local(__file__).dirpath(
        "..", "docs", "quick_start", "basak_tabulate_lammps")
    srcdir.join("UO2.lmpstruct").copy(tmpdir.join("UO2.lmpstruct"))

    if aspotfile.basename == "basak_table_form.aspot":
        input_file = py.path.local(_getLAMMPSResourceDirectory()).join(
            "basak_energy_table_form.lmpin")
        expect_e = pytest.approx(-172.839, abs=1e-3)
    else:
        input_file = py.path.local(
            _getLAMMPSResourceDirectory()).join("basak_energy.lmpin")
        expect_e = pytest.approx(-172.924, abs=1e-3)

    input_file.copy(tmpdir.join("calc_energy.lmpin"))

    # Generate table file
    config = atsim.potentials.config.Configuration()
    tabulation = config.read(aspotfile.open())

    with tmpdir.join("Basak.lmptab").open("w") as outfile:
        tabulation.write(outfile)

    runLAMMPS(cwd=tmpdir.strpath)
    actual_e = extractLAMMPSEnergy(cwd=tmpdir.strpath)

    assert expect_e == actual_e


morelon_files = py.path.local(__file__).dirpath(
    "..", "docs", "user_guide", "example_files").listdir("morelon*.aspot")


@needsGULP
@pytest.mark.parametrize("charges", [[-1.613626, 3.227252]])
@pytest.mark.parametrize("aspot", morelon_files)
def test_morelon_files(aspot, gulp_uo2_energy_fixture):
    tmpdir = gulp_uo2_energy_fixture
    # aspot = py.path.local(__file__).dirpath("..", "docs", "user_guide", "example_files").join("morelon.aspot")

    CPT = atsim.potentials.config.ConfigParserOverrideTuple
    # overrides = [ CPT("Tabulation", "target", "GULP"), CPT("Tabulation", "cutoff", "10.0"), CPT("Tabulation", "nr", "1001")]
    overrides = [CPT("Tabulation", "target", "GULP")]

    # import io
    with aspot.open('r', encoding='utf-8') as infile:
        cp = atsim.potentials.config.ConfigParser(infile, overrides=overrides)
        tabulation = atsim.potentials.config.Configuration().read_from_parser(cp)

    with tmpdir.join("potentials.lib").open("w", encoding='utf-8') as tabfile:
        tabulation.write(tabfile)

    expect = pytest.approx(-263.60598484)
    actual_energy = gulp_uo2_energy_fixture.energy()

    assert expect == actual_energy


@needsLAMMPS
@pytest.mark.parametrize(("evaluate_lmpin", "aspot_filename", "tab_filename", "expect_e", "a_flag", "b_flag"),
                         [
    ("standard_evaluate.lmpin", "standard_eam.aspot",
     "standard_eam.eam", 24.00, True, False),
    ("standard_evaluate.lmpin", "standard_eam.aspot",
        "standard_eam.eam", 131.8822, False, True),
    ("standard_evaluate.lmpin", "standard_eam.aspot",
        "standard_eam.eam", 24.00 + 131.8822, True, True),

    ("finnis_sinclair_evaluate.lmpin", "finnis_sinclair_eam.aspot",
        "finnis_sinclair.eam.fs", 24.00, True, False),
    ("finnis_sinclair_evaluate.lmpin", "finnis_sinclair_eam.aspot",
        "finnis_sinclair.eam.fs", 209.137, False, True),
    ("finnis_sinclair_evaluate.lmpin", "finnis_sinclair_eam.aspot",
        "finnis_sinclair.eam.fs", 24.00 + 209.137, True, True)
])
def test_user_guide_eam(tmpdir, evaluate_lmpin, aspot_filename, tab_filename, expect_e, a_flag, b_flag):

    # Copy files in from example directory
    srcdir = py.path.local(__file__).dirpath(
        "..", "docs", "user_guide", "example_files")
    srcdir.join("toy_structure.lmpstruct").copy(
        tmpdir.join("toy_structure.lmpstruct"))

    input_file = srcdir.join(evaluate_lmpin)
    input_file.copy(tmpdir.join("calc_energy.lmpin"))

    # Generate table file
    config = atsim.potentials.config.Configuration()

    aspotfile = srcdir.join(aspot_filename)

    a_embed = 'as.zero'
    b_embed = 'as.zero'

    if a_flag:
        a_embed = 'as.polynomial 0 1'

    if b_flag:
        b_embed = 'as.polynomial 0 1'

    overrides = [
        atsim.potentials.config.ConfigParserOverrideTuple(
            u'EAM-Embed', 'A', a_embed),
        atsim.potentials.config.ConfigParserOverrideTuple(
            u'EAM-Embed', 'B', b_embed)
    ]

    cp = atsim.potentials.config.ConfigParser(aspotfile.open(), overrides)
    tabulation = config.read_from_parser(cp)

    # tabulation = config.read(aspotfile.open())

    with tmpdir.join(tab_filename).open("w") as outfile:
        tabulation.write(outfile)

    expect_e = pytest.approx(expect_e, abs=1e-3)
    runLAMMPS(cwd=tmpdir.strpath)
    actual_e = extractLAMMPSEnergy(cwd=tmpdir.strpath)

    assert expect_e == actual_e
