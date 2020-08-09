import pytest

import py.path

import distutils.spawn
GULP_FOUND = distutils.spawn.find_executable("gulp-5.0")

needsGULP = pytest.mark.skipif(not GULP_FOUND, reason = "GULP binary not found")

def _getResourceDirectory():
  """Returns path to resources used by this test module (currently assumed to be sub-directory
  of test module called resources)"""
  return py.path.local(__file__).dirpath('gulp_resources')

@pytest.fixture
def charges():
  return [-2, 4]

@pytest.fixture
def gulp_uo2_energy_fixture(tmpdir, charges):
  from io import StringIO

  def run(charges):
    resource_path = py.path.local(_getResourceDirectory())
    gulpin = resource_path.join("uo2.gin")
    template = gulpin.open(encoding = 'utf-8').read()

    ocharge, ucharge = charges
    charges = {"U_charge" : ucharge, "O_charge" : ocharge}
    
    outfile = StringIO()
    infile = StringIO(template.format(**charges))
    infile.seek(0)

    runGULP(infile, outfile, cwd = tmpdir.strpath)
    outfile.seek(0)

    return outfile
  
  def energy():
    outfile = run(charges)
    energy = extractGULPEnergy(outfile)
    return energy

  tmpdir.energy = energy

  return tmpdir


def runGULP(infile, outfile, cwd = None):  
  import subprocess
  popen = subprocess.Popen("gulp-5.0", cwd = cwd, shell = True, stdin = subprocess.PIPE, stdout = subprocess.PIPE, close_fds=True)
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

