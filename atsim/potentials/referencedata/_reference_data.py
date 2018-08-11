from ._data import reference_data, Element_Data


class Reference_Data_Exception(Exception):
  pass

class Unknown_Species_Exception(Reference_Data_Exception):
  pass

class Unknown_Property_Exception(Reference_Data_Exception):
  pass

class Reference_Data(object):
  """Class providing data about atomic species"""
  

  def __init__(self, extra_data = {}):
    """Create object for looking up reference data about atomic species.
    
    :param extra_data: Override defaults or provide additional data. This argument is a dictionary of dictionaries:
      `{SPECIES_A : SPECIES_DATA}`

        Where `SPECIES_DATA` maps property names to property values, e.g. `{ 'atomic_mass` : 1.234}`.
    """
    self.extra_data = extra_data

  def get(self, species, property_name):
    """Get a property value for a given species.

    :param species: Species label.
    :param property_name: Propety identifier.

    :returns: Property value for given combination of species and property name."""
    species_dat = reference_data.get(species, None)

    if species_dat is None and not species in self.extra_data:
      raise Unknown_Species_Exception(species)
    elif species_dat is None and species in self.extra_data:
      species_dat = dict(self.extra_data[species])
    else:
      species_dat = species_dat._asdict()
      species_dat.update(self.extra_data.get(species, {}))

    if not property_name in species_dat:
      raise Unknown_Property_Exception("Property '{}' not found for species '{}'".format(property_name, species))
    return species_dat[property_name]