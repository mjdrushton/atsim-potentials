# from future import standard_library
# standard_library.install_aliases()
import unittest
import py.path
import pytest

import os
import distutils.spawn
LAMMPS_FOUND = distutils.spawn.find_executable('lammps')

# needsLAMMPS = unittest.skipIf(not LAMMPS_FOUND, "LAMMPS not available")
needsLAMMPS = pytest.mark.skipif(not LAMMPS_FOUND, reason="LAMMPS not available")

def _getResourceDirectory():
  """Returns path to resources used by this test module (currently assumed to be sub-directory
  of test module called resources)"""
  return os.path.join(os.path.dirname(__file__), 'lammps_resources')

@pytest.fixture
def lammps_run_fixture(tmpdir):
  resource_path = py.path.local(_getResourceDirectory())
  lmpin = resource_path.join("calc_energy.lmpin")
  
  import shutil
  shutil.copyfile(lmpin.strpath, tmpdir.join("calc_energy.lmpin").strpath)
  shutil.copyfile(resource_path.join("pair_2angs.lmpstruct").strpath,
    tmpdir.join("structure.lmpstruct").strpath)

  return tmpdir

@pytest.fixture
def lammps_run_fluorite_fixture(tmpdir):
  resource_path = py.path.local(_getResourceDirectory())
  lmpin = resource_path.join("calc_energy.lmpin")
  
  import shutil
  shutil.copyfile(lmpin.strpath, tmpdir.join("calc_energy.lmpin").strpath)
  shutil.copyfile(resource_path.join("CeO2-single_cell.lmpstruct").strpath,
    tmpdir.join("structure.lmpstruct").strpath)

  return tmpdir

def extractLAMMPSEnergy(cwd = None):
  if cwd:
    outpath = os.path.join(cwd, "out.lmpout")
  else:
    outpath = "out.lmpout"

  with open(outpath) as infile:
    for line in infile:
      line = line[:-1]
      if line.startswith('ENERGY:'):
        tokens = line.split(':')
        energy = float(tokens[1])
        return energy

def runLAMMPS(cwd = None):
  import subprocess
  output = subprocess.check_output("lammps -in calc_energy.lmpin -log out.lmpout -echo none", shell = True, cwd = cwd)