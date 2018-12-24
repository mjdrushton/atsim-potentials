import pytest

import distutils
GULP_FOUND = distutils.spawn.find_executable("gulp-5.0")

needsGULP = pytest.mark.skipif(not GULP_FOUND, reason = "GULP binary not found")

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

