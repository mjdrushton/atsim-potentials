import logging

from .._eam_potential import EAMPotential
from ._common import ConfigurationException
from ._potential_form_builder import Potential_Form_Builder

from ..referencedata import reference_data

class EAM_Potential_Builder(object):
  """Uses the output ConfigParser.eam_density and .eam_embed propoerties and builds
  EAMPotential instances"""


  def __init__(self, cp, potential_form_registry, modifier_registry):
    """:param cp: atsim.potentials.config.ConfigParser instance.
    :param potential_form_register: Potential_Form_Registry
    :param modifier_register: Modifier_Registry"""
    self._potlist = self._init_eampotentials(cp, potential_form_registry, modifier_registry)

  def _init_eampotentials(self, cp, potential_form_registry, modifier_registry):
    pots = []

    # Gather the embedding functions
    embed = cp.eam_embed
    # Gather the density functions
    density = cp.eam_density

    # Check that we have the same set of species in embed and density
    embed_species = set([row.species for row in embed])
    density_species = set([row.species for row in density])

    diff = embed_species ^ density_species

    if diff:
      embed_species = ",".join(sorted(embed_species))
      density_species = ",".join(sorted(density_species))
      diff_species = ",".join(sorted(diff))
      errmsg = "During EAM tabulation species defined for density function do not match those for embedding functions. Density species: {density}. Embed species: {embed}. Difference: {difference}"
      errmsg = errmsg.format(density = density_species,
        embed_species = embed_species,
        difference = diff_species)
      raise ConfigurationException(errmsg)

    # Convert the embed and density tuple lists into
    # maps relating species to a potential function.
    potential_form_builder = Potential_Form_Builder(potential_form_registry, modifier_registry)
    embed_dict = self._to_potential_form_dict(embed, potential_form_builder)
    density_dict = self._to_potential_form_dict(density, potential_form_builder)

    # Now instantiate EAM potential objects
    potlist = []
    for species in embed_dict:
      pot = self._create_eam_potential(species, embed_dict, density_dict)
      potlist.append(pot)
    return potlist

  def _to_potential_form_dict(self, tuple_list, potential_form_builder):
    d = {}

    for t in tuple_list:
      species = t.species
      #Create potential function
      pot_func = potential_form_builder.create_potential_function(t.potential_form_instance)
      d[species] = pot_func
    return d

  def _get_mass(self, species):
    try:
      return reference_data[species].atomic_mass
    except KeyError:
      raise ConfigurationException("Could not find atomic mass for species: {}".format(species))

  def _get_atomic_number(self, species):
    try:
      return reference_data[species].atomic_number
    except KeyError:
      raise ConfigurationException("Could not find atomic number for species: {}".format(species))

  def _get_lattice_constant(self, species):
    return 0.0
    
  def _get_lattice_type(self, species):
    return 'fcc'

  def _create_eam_potential(self, species, embed_dict, density_dict):
    embedding_function = embed_dict[species]
    density_function = density_dict[species]

    atomic_number = self._get_atomic_number(species)
    mass = self._get_mass(species)
    lattice_constant = self._get_lattice_constant(species)
    lattice_type = self._get_lattice_type(species)

    eam_pot = EAMPotential(species, 
      atomic_number, 
      mass, 
      embedding_function, density_function,
      lattice_constant,
      lattice_type)
    
    return eam_pot

  @property
  def eam_potentials(self):
    return self._potlist