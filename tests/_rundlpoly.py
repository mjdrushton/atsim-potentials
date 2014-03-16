import unittest

import distutils
DLPOLY_FOUND = distutils.spawn.find_executable('DLPOLY.Z')

needsDLPOLY = unittest.skipIf(not DLPOLY_FOUND, "DLPOLY not available")


def runDLPoly():
  import commands
  commands.getoutput('DLPOLY.Z')

def extractDLPOLYEnergy():
  with open('STATIS') as infile:
    # import pdb;pdb.set_trace()
    infile.next()
    infile.next()
    infile.next()
    line = infile.next()
    tokens = line.split()
    engcfg = float(tokens[2])
  return engcfg
