#! /usr/bin/env python

import atsim_potentials
from atsim_potentials import potentialforms

def main():
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
    morse_OU = potentialforms.morse(0.577190, 1.6500, 2.36900)

    # Compose the buckingham and morse functions into a single function
    # using the atsim_potentials.plus() function
    f_OU = atsim_potentials.plus(buck_OU, morse_OU)

    # Construct list of Potential objects
    potential_objects = [
        atsim_potentials.Potential('O', 'O', f_OO),
        atsim_potentials.Potential('U', 'U', f_UU),
        atsim_potentials.Potential('O', 'U', f_OU)
    ]

    # Tabulate into file called TABLE
    # using short-range cutoff of 6.5 Angs with grid
    # increment of 1e-3 Angs (6500 grid points)
    with open('TABLE', 'wb') as outfile:
        atsim_potentials.writePotentials(
           'DL_POLY',
           potential_objects,
           6.5, 6500,
           out = outfile)

if __name__ == '__main__':
    main()
