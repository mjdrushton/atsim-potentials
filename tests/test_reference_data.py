import pytest

from atsim.potentials.referencedata import Reference_Data, Unknown_Species_Exception, Unknown_Property_Exception


def test_mass():
  rd = Reference_Data()
  actual = rd.get("Gd", "atomic_mass")
  assert pytest.approx(157.25) == actual

def test_mass_override():
  rd = Reference_Data({"Gd" : {"atomic_mass" : 924.0 }})
  actual = rd.get("Gd", "atomic_mass")
  assert pytest.approx(924) == actual
  assert pytest.approx(195.078) == rd.get("Pt", "atomic_mass")

def test_unknown_species():
  rd = Reference_Data()

  with pytest.raises(Unknown_Species_Exception):
    rd.get("Bl", 'atomic_mass')

def test_unknown_species():
  rd = Reference_Data()

  with pytest.raises(Unknown_Property_Exception):
    rd.get("Ba", 'atomic_massive')