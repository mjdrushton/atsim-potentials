import io
import sys

from atsim.potentials.config import (ConfigParser, ConfigParserOverrideTuple,
                                     Configuration)

potable_input = """
[Tabulation]
target : GULP
cutoff : 10.0
dr : 0.01

[Pair]
Si-O : spline(
                    as.zbl 14 8 
              >=0.8 
                    exp_spline 
              >=1.4 
                    as.buck 18003.7572 0.205204 133.5381 )
"""


def main():
    # Make a file like object from the potable input string given above.
    potable_input_file = io.StringIO(potable_input)

    # Create a Configuration() object and read input from the input file.
    configuration = Configuration()

    # This example shows how to override and add items to potable input before it
    # is passed to the Configuration object.
    #    The tabulation target will be change to 'LAMMPS'
    #    The cutoff will be reduced to 6.5
    #    An additional pair-interaction will be given for O-O

    cp = ConfigParser(potable_input_file,
                      overrides=[
                          ConfigParserOverrideTuple(
                              "Tabulation", "target", "LAMMPS"),
                          ConfigParserOverrideTuple(
                              "Tabulation", "cutoff", "6.5")
                      ],
                      additional=[
                          ConfigParserOverrideTuple(
                              "Pair", "O-O", "as.buck 444.7686 0.402 0.0")
                      ])

    # Create the tabulation by passing the Config_Parser object to the Configuration.read_from_parser method.
    tabulation = configuration.read_from_parser(cp)

    # Now write tabulation to console
    tabulation.write(sys.stderr)


if __name__ == "__main__":
    main()
