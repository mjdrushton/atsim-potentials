#! /usr/bin/env python
import math

from atsim.potentials import EAMPotential, Potential
from atsim.potentials.eam_tabulation import SetFL_EAMTabulation


def embed(rho):
    return -math.sqrt(rho)


def density(rij):
    if rij == 0:
        return 0.0
    return (2.928323832 / rij) ** 6.0


def pair_AgAg(rij):
    if rij == 0:
        return 0.0
    return (2.485883762/rij) ** 12


def main():
    # Create EAMPotential
    eam_potentials = [EAMPotential("Ag", 47, 107.8682, embed, density)]
    pair_potentials = [Potential('Ag', 'Ag', pair_AgAg)]

    cutoff_rho = 50.0
    nrho = 50000

    cutoff = 12.0
    nr = 12000

    tabulation = SetFL_EAMTabulation(
        pair_potentials,
        eam_potentials,
        cutoff, nr,
        cutoff_rho, nrho)

    with open("Ag.eam.alloy", 'w') as outfile:
        tabulation.write(outfile)


if __name__ == "__main__":
    main()
