#! /usr/bin/env python

from atsim.potentials import potentialforms
import atsim.potentials

def main():

    bks_buck = potentialforms.buck(18003.7572, 1.0/4.87318, 133.5381)
    bks_coul = potentialforms.coul(2.4, -1.2)

    bks = atsim.potentials.plus(bks_buck, bks_coul)

    zbl = potentialforms.zbl(14, 8)


    spline = atsim.potentials.SplinePotential( zbl, bks_buck, 0.8, 1.4)

    atsim.potentials.plot( 'bks_buck.dat', 0.1, 10.0, bks_buck, 5000)
    atsim.potentials.plot( 'bks_coul.dat', 0.1, 10.0, bks_coul, 5000)
    atsim.potentials.plot( 'bks.dat', 0.1, 10.0, bks, 5000)
    atsim.potentials.plot( 'zbl.dat', 0.1, 10.0, zbl, 5000)
    atsim.potentials.plot( 'spline.dat', 0.1, 10.0, spline, 5000)

    bks_SiO = atsim.potentials.Potential('Si', 'O', spline)
    with open('bks_SiO.lmptab', 'wb') as outfile:
        atsim.potentials.writePotentials('LAMMPS', [bks_SiO], 10.0, 5000, out = outfile)


if __name__ == '__main__':
    main()
