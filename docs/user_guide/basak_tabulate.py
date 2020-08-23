#! /usr/bin/env python

from atsim.potentials import Potential, plus, potentialforms
from atsim.potentials.pair_tabulation import DLPoly_PairTabulation


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
    f_OU = plus(buck_OU, morse_OU)

    # Construct list of Potential objects
    potential_objects = [
        Potential('O', 'O', f_OO),
        Potential('U', 'U', f_UU),
        Potential('O', 'U', f_OU)
    ]
    return potential_objects


def main():
    potential_objects = makePotentialObjects()
    # Tabulate into file called TABLE
    # using short-range cutoff of 6.5 Angs with grid
    # increment of 1e-3 Angs (6500 grid points)

    tabulation = DLPoly_PairTabulation(potential_objects,
                                       6.5, 6500)

    with open('TABLE', 'w') as outfile:
        tabulation.write(outfile)


if __name__ == '__main__':
    main()
