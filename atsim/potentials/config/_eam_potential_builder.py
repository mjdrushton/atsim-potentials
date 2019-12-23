import logging

from .._eam_potential import EAMPotential
from ._common import ConfigurationException
from ._potential_form_builder import Potential_Form_Builder
from ..potentialforms import zero

from ..referencedata import Reference_Data, Reference_Data_Exception

class EAM_Potential_Builder(object):
  """Uses the output ConfigParser.eam_density and .eam_embed properties and builds
  EAMPotential instances"""

  def __init__(self, cp, potential_form_registry, modifier_registry, reference_data = Reference_Data(), add_undefined = True):
    """:param cp: atsim.potentials.config.ConfigParser instance.
    :param potential_form_register: Potential_Form_Registry
    :param modifier_register: Modifier_Registry
    :param reference_data: Reference_Data object used to provide atomic number, masses and default lattices.
    :param add_undefined: Add null embedding functions and density function for missing interactions. This allows underspecified systems to be tabulatd more easily."""
    self._reference_data = reference_data
    self.add_undefined = add_undefined
    self._potlist = self._init_eampotentials(cp, potential_form_registry, modifier_registry)

  def _init_eampotentials(self, cp, potential_form_registry, modifier_registry):
    # Gather the embedding functions
    embed = self._extract_embed(cp)
    # Gather the density functions
    density = self._extract_density(cp)

    # Check that we have the same set of species in embed and density
    embed_species = self._embed_species(embed) 
    density_species = self._density_species(density)

    diff = embed_species ^ density_species

    if diff and not self.add_undefined:
      embed_species = ",".join(sorted(embed_species))
      density_species = ",".join(sorted(density_species))
      diff_species = ",".join(sorted(diff))
      errmsg = "During EAM tabulation species defined for density function do not match those for embedding functions. Density species: {density}. Embed species: {embed_species}. Difference: {difference}"
      errmsg = errmsg.format(density = density_species,
        embed_species = embed_species,
        difference = diff_species)
      raise ConfigurationException(errmsg)

    # Convert the embed and density tuple lists into
    # maps relating species to a potential function.
    potential_form_builder = Potential_Form_Builder(potential_form_registry, modifier_registry)
    embed_dict = self._embed_to_potential_form_dict(embed, potential_form_builder)
    density_dict = self._density_to_potential_form_dict(density, potential_form_builder)

    if self.add_undefined:
      self._add_null_functions(cp, embed_dict, density_dict)

    # Now instantiate EAM potential objects
    potlist = []
    for species in embed_dict:
      pot = self._create_eam_potential(species, embed_dict, density_dict)
      potlist.append(pot)
    return potlist

  def _add_null_functions(self, cp, embed_dict, density_dict):
    # Add default no-op embedding and density functions for species defined in configuration file but that
    # have not been included in embed_dict and density_dict at this point in potential building
    self._add_null_embedding_functions(cp, embed_dict, density_dict)
    self._add_null_density_functions(cp, embed_dict, density_dict)

  def _pp_species(self, cp):
    pair_species = set()
    for pp_tup in cp.pair:
      pair_species.add(pp_tup.species.species_a)
      pair_species.add(pp_tup.species.species_b)
    return pair_species

  def _add_null_embedding_functions(self, cp, embed_dict, density_dict):
    # Get set of currently defined embedding functions
    defined = set(embed_dict.keys())

    # Get set of species defined for pair potentials
    # pair_species = self._pp_species(cp)
    
    # Get set of species defined by density functions
    density = self._extract_density(cp)
    density_species = self._density_species(density)

    # null_embed_species = (pair_species | density_species) - defined
    null_embed_species = density_species - defined

    # Create the zero functions for null_embed_species.
    null = zero()
    for s in null_embed_species:
      embed_dict[s] = null


  def _add_null_density_functions(self, cp, embed_dict, density_dict):
    embed_species = set(embed_dict.keys())
    # pair_species = self._pp_species(cp)
    density = self._extract_density(cp)
    density_species = self._density_species(density)
    # all_species = embed_species | pair_species | density_species
    all_species = embed_species | density_species

    null = zero()
    for s in all_species:
      other_dict = density_dict.setdefault(s, null)


  def _extract_embed(self, cp):
    return cp.eam_embed

  def _extract_density(self, cp):
    return cp.eam_density

  def _embed_species(self, embed):
    return set([row.species for row in embed])

  def _density_species(self, density):
    return set([row.species for row in density])

  def _embed_to_potential_form_dict(self, tuple_list, potential_form_builder):
    return self._to_potential_form_dict(tuple_list, potential_form_builder)

  def _density_to_potential_form_dict(self, tuple_list, potential_form_builder):
    return self._to_potential_form_dict(tuple_list, potential_form_builder)

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
      return self._reference_data.get(species, 'atomic_mass')
    except Reference_Data_Exception:
      raise ConfigurationException("Could not find atomic mass for species: {}".format(species))

  def _get_atomic_number(self, species):
    try:
      return self._reference_data.get(species, 'atomic_number')
    except Reference_Data_Exception:
      raise ConfigurationException("Could not find atomic number for species: {}".format(species))

  def _get_lattice_constant(self, species):
    try:
      return self._reference_data.get(species, 'lattice_constant')
    except Reference_Data_Exception:
      return 0.0
    
  def _get_lattice_type(self, species):
    try:
      return self._reference_data.get(species, 'lattice_type')
    except Reference_Data_Exception:
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


class EAM_Potential_Builder_FS(EAM_Potential_Builder):
  """EAM_Potential_Builder for Finnis-Sinclair style potentials.

  The major difference between FS potentials and standard potentials
  is that density functions are defined in EAMPotential as  dictionaries
  as densities are specific to interacting pairs of species"""

  def _extract_density(self, cp):
    return cp.eam_density_fs

  def _density_species(self, density):
    species_list = []
    for row in density:
      species_list.append(row.species.from_species)
      species_list.append(row.species.to_species)
    return set(species_list)

  def _density_to_potential_form_dict(self, density, potential_form_builder):
    outdict = {}
    for d in density:
      f_species = d.species.from_species
      t_species = d.species.to_species
      
      pot_func = potential_form_builder.create_potential_function(d.potential_form_instance)

      add_to = outdict.setdefault(f_species, {})
      
      if t_species in add_to:
        raise ConfigurationException("Duplicate density function found for {}".format(d.species))

      add_to[t_species] = pot_func
    return outdict

  def _add_null_density_functions(self, cp, embed_dict, density_dict):
    embed_species = set(embed_dict.keys())
    # pair_species = self._pp_species(cp)
    density = self._extract_density(cp)
    density_species = self._density_species(density)
    # all_species = embed_species | pair_species | density_species
    all_species = embed_species | density_species

    null = zero()
    for s in all_species:
      other_dict = density_dict.setdefault(s, {})
      for o in all_species:
        other_dict.setdefault(o, null)