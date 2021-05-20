class EAMPotential(object):
  """Class used to describe a particular species within EAM potential models.

  This class is a container for the functions and attributes necesary for describing the many-body component of an
  Embedded Atom potential Model.


  """

  def __init__(self,
              species,
              atomicNumber,
              mass,
              embeddingFunction,
              electronDensityFunction,
              latticeConstant = 0.0,
              latticeType = 'fcc'):
    """Create EAMPotential object.

    **Note on electronDensityFunction parameter**

    In the conventional Embedded Atom Method each ``EAMPotential`` encapsulates a single density function
    (density functions are mixed by the simulation code for alloy systems). In such cases the ``electronDensityFunction``
    constructor parameter is a single python callable.

    The Finnis-Sinclair extension of the Embdedded Atom Method allows density functions for specific pairs of species
    to be used. This requires a callable to be given for each species pair. This is achieved by passing a dictionary
    of callables/functions to the ``electronDensityFunction`` constructor parameter. This has the general form:

    .. code::

        { SPECIES_1 : DENSITY_FUNCTION_1,
          SPECIES_2 : DENSITY_FUNCTION_2,
          ...
          SPECIES_N : DENSITY_FUNCTION_N}


    Where the keys ``SPECIES_N`` represent the species of atoms surrounding a central atom with the ``species``
    passed to the ``EAMPotential`` constructor and ``DENSITY_FUNCTION_N`` is a unary function that takes the separation
    between the central atom and the surrounding atom and returns the electron density due to the surrounding atom.

    **Example:**
    For a binary system containing Al and Cu, the ``EAMPotential`` constructor may be invoked as follows to create an
    object representing Al for use with
    the Finnis-Sinclair style EAM tabulation functions:

    .. code::

      alPotential = EAMPotential(
         "Al",
         13,
         26.98,
         embeddingFunction,
         {"Al" : density_Al_Al, # specifies Al-Al density function
          "Cu" : density_Al_Cu, # specifies Al-Cu density function
         })



    .. seealso::

      * For working examples of how to use ``EAMPotential`` with the conventional Embedded Atom Model:

        * :ref:`eam_example_1`
        * :ref:`eam_example_2a`
        * :ref:`eam_example_2b`

      * For complete examples using the Finnis-Sinclair version of the method:

        * :ref:`eam_example_3a`
        * :ref:`eam_example_3b`


    :param species: Species label
    :param atomicNumber: Atomic number
    :param mass: Atomic mass
    :param embeddingFunction: Python callable that retuns embedding values as function of density
    :param electronDensityFunction: Python callable that returns density as function of separation in angstroms.
       For write functions that accept multiple density functions (e.g. :func:`~atsim.potentials.writeSetFLFinnisSinclair` and
       :func:`~atsim.potentials.writeTABEAMFinnisSinclair`) a dictionary mapping species to a function object is expected for this
       parameter. See above.
    :param latticeConstant: Lattice constant for this species
    :param latticeType: Lattice type e.g. 'fcc' for this"""

    #: Atomic species represented by this object.
    self.species = species

    #: Atomic number (unused by most simulation codes)
    self.atomicNumber = atomicNumber

    #: Atomic mass (unused by most simulation codes)
    self.mass = mass

    #: Python callable taking single argument (density) and returns energy assumed to be in eV.
    self.embeddingFunction = embeddingFunction

    #: Python callable that returns density as function of separation in angstroms.
    #  For write functions that accept multiple density functions (e.g. writeSetFLFinnisSinclair() )
    #  a dictionary mapping element to a function object is expected for this parameter.
    #  See :meth:`atsim.potentials.EAMPotential.__init__` for more.
    self.electronDensityFunction = electronDensityFunction

    #: Lattice constant for species (unused by most simulation codes).
    self.latticeConstant = latticeConstant

    #: Lattice type e.g. FCC (unused by most simulation codes).
    self.latticeType = latticeType

  def embeddingValue(self, density):
    """Method that returns energy for given electron density.

    This method simply passes ``density`` to the callable stored in the ``embeddingFunction`` and returns its value.

    :param density: Electron density.
    :type density: float
    :return: Energy for given density (as given by ``self.embeddingFunction``).
    :rtype: float"""
    return self.embeddingFunction(density)

  def electronDensity(self, separation):
    """Gives the 'electron' density for an atom separated from current species by ``separation``.

    This is a pass-through method to callable stored in current instance's ``electronDensityFunction`` attribute.

    .. :warning:: In cases where the ``electronDensityFunction`` attribute has been set to a dictionary, (for instance, when density functions
      specific to particular pairs of species are required), this method should not be used

    :param separation: Separation (in angstroms) between atom represented by this object and another atom.
    :type separation: float.

    :return: Contribution to electron density due to given pair separation.
    :rtype: float."""
    return self.electronDensityFunction(separation)
