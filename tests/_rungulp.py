import pytest

import pathlib
import types

import distutils.spawn
GULP_FOUND = distutils.spawn.find_executable("gulp-5.2")

needsGULP = pytest.mark.skipif(not GULP_FOUND, reason = "GULP binary not found")

def _getResourceDirectory():
  """Returns path to resources used by this test module (currently assumed to be sub-directory
  of test module called resources)"""
  return pathlib.Path(__file__).parent / 'gulp_resources'

@pytest.fixture
def charges():
  return [-2, 4]

@pytest.fixture
def gulp_uo2_energy_fixture(tmp_path, charges):
  from io import StringIO

  class Gulp_Directory(pathlib.Path):

    def __new__(cls, p: pathlib.Path, charges):
      self = pathlib.Path.__new__(cls, p)
      self.charges = charges
      return self

    def run(self, charges):
      resource_path = pathlib.Path(_getResourceDirectory())
      gulpin = resource_path / "uo2.gin"
      template = gulpin.open(encoding = 'utf-8').read()

      ocharge, ucharge = charges
      charges = {"U_charge" : ucharge, "O_charge" : ocharge}
      
      outfile = StringIO()
      infile = StringIO(template.format(**charges))
      infile.seek(0)

      runGULP(infile, outfile, cwd = tmp_path)
      outfile.seek(0)

      return outfile
    
    def energy(self):
      outfile = self.run(self.charges)
      energy = extractGULPEnergy(outfile)
      return energy

  ret_path = Gulp_Directory(tmp_path, charges)
  
  return ret_path


def runGULP(infile, outfile, cwd = None):  
  import subprocess
  popen = subprocess.Popen("gulp-5.2", cwd = cwd, shell = True, stdin = subprocess.PIPE, stdout = subprocess.PIPE, close_fds=True)
  stdout, stderr = popen.communicate(infile.read().encode("utf-8"))
  outfile.write(stdout.decode())

def extractGULPEnergy(infile):
  markers= ['*  Output for configuration   1',
  '  Components of energy :',
  '  Total lattice energy       =']

  for find_line in markers:
    found = False
    for line in infile:
      if line.startswith(find_line):
        found = True
        break
    if not found:
      pytest.fail("Could not extract GULP energy")

  line = line.strip()
  assert line.endswith("eV")

  line = line[:-2]
  line, E = line.split("=")

  E = float(E)
  return E

