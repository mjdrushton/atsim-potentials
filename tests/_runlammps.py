import unittest

import distutils
LAMMPS_FOUND = distutils.spawn.find_executable('lammps')

needsLAMMPS = unittest.skipIf(not LAMMPS_FOUND, "LAMMPS not available")

def extractLAMMPSEnergy():
  with open('out.lmpout') as infile:
    for line in infile:
      line = line[:-1]
      if line.startswith('ENERGY:'):
        tokens = line.split(':')
        energy = float(tokens[1])
        return energy

def runLAMMPS():
  import commands
  output = commands.getoutput("lammps -in calc_energy.lmpin -log out.lmpout -echo none")

