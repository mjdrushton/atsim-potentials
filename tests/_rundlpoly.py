# from future import standard_library
# standard_library.install_aliases()
from builtins import next
import unittest
import os

import distutils
DLPOLY_FOUND = distutils.spawn.find_executable('DLPOLY.Z')

needsDLPOLY = unittest.skipIf(not DLPOLY_FOUND, "DLPOLY not available")


def runDLPoly(cwd = None):
  import subprocess
  subprocess.check_output('DLPOLY.Z', shell = True, cwd = cwd)

def extractDLPOLYEnergy(cwd = None):
  if cwd:
    statispath = os.path.join(cwd, "STATIS")
  else:
    statispath = "STATIS"

  with open(statispath) as infile:
    next(infile)
    next(infile)
    next(infile)
    line = next(infile)
    tokens = line.split()
    engcfg = float(tokens[2])
  return engcfg
