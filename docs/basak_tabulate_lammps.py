#! /usr/bin/env python

import atsim.potentials
from atsim.potentials import potentialforms

def makePotentialObjects():
    # O-O Interaction:
    # Buckingham
    # A = 1633.00510, rho = 0.327022, C = 3.948790
    f_OO = potentialforms.buck(1633.00510, 0.327022, 3.948790)

    # U-U Interaction:
    # Buckingham
    # A = 294.640000, rho = 0.327022, C = 0.0
    f_UU = potentialforms.buck(294.640000, 0.327022, 0.0)

    # O-U Interaction
    # Buckingham + Morse.
    # Buckingham:
    # A = 693.648700, rho = 693.648700, C = 0.0
    # Morse:
    # D0 = 0.577190, alpha = 1.6500, r0 = 2.36900
    buck_OU = potentialforms.buck(693.648700, 0.327022, 0.0)
    morse_OU = potentialforms.morse(1.6500, 2.36900, 0.577190)


    # Compose the buckingham and morse functions into a single function
    # using the atsim.potentials.plus() function
    f_OU = atsim.potentials.plus(buck_OU, morse_OU)

    # Construct list of Potential objects
    potential_objects = [
        atsim.potentials.Potential('O', 'O', f_OO),
        atsim.potentials.Potential('U', 'U', f_UU),
        atsim.potentials.Potential('O', 'U', f_OU)
    ]
    return potential_objects

def main():
    potential_objects = makePotentialObjects()
    # Tabulate into file called Basak.lmptab
    # using short-range cutoff of 6.5 Angs with grid
    # increment of 1e-3 Angs (6500 grid points)
    with open('Basak.lmptab', 'wb') as outfile: # <-- Filename changed from 'TABLE'
        atsim.potentials.writePotentials(
           'LAMMPS', # <-- This has been changed from 'DL_POLY'
           potential_objects,
           6.5, 6500,
           out = outfile)

if __name__ == '__main__':
    main()
