from future import standard_library
standard_library.install_aliases()
from builtins import next
import unittest

import distutils
DLPOLY_FOUND = distutils.spawn.find_executable('DLPOLY.Z')

needsDLPOLY = unittest.skipIf(not DLPOLY_FOUND, "DLPOLY not available")


def runDLPoly():
  import subprocess
  subprocess.getoutput('DLPOLY.Z')

def extractDLPOLYEnergy():
  with open('STATIS') as infile:
    # import pdb;pdb.set_trace()
    next(infile)
    next(infile)
    next(infile)
    line = next(infile)
    tokens = line.split()
    engcfg = float(tokens[2])
  return engcfg
