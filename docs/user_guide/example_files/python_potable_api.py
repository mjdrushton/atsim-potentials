import io

from atsim.potentials.config import Configuration

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
    # ... Configuration is a factory class for PairTabulation and EAMTabulation
    #     objects. In the current case it will return a GULP_PairTabulation object.
    tabulation = configuration.read(potable_input_file)

    # The potable input defines a single pair potential.
    # Potential objects are accessible from the tabulation object through
    # its .potentials property.
    potential_Si_O = tabulation.potentials[0]

    # The potential-form for this interaction is now accessed.
    multirange_potentialform = potential_Si_O.potentialFunction

    # The potential-forms created from potable input are Multi_Range_Potential_Form
    # objects. This is true even if only one range is defined, as is the case here.
    #
    # Let's get hold of the spline potential form through the Multi_Range_Potential_Form
    # .range_defns property (which returns a list of MultiRangeDefinitionTuple)
    #
    spline_potentialform = multirange_potentialform.range_defns[0].potential_form

    # Now let's get hold of the spline coefficients
    spline_coefficients = spline_potentialform.splineCoefficients

    print("Spline coefficients are: {}".format(spline_coefficients))


if __name__ == "__main__":
    main()
