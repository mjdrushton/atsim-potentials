import pytest

from ._common import _get_dlpoly_resource_dir, _get_lammps_resource_dir

from atsim.potentials.tools.potable import _query_actions
from atsim.potentials.config import ConfigParser


def test_list_item_labels():
  expect = [
    "Tabulation:target",
    "Tabulation:nr",
    "Tabulation:dr",
    "Tabulation:nrho",
    "Tabulation:drho",
    "Potential-Form:buck_morse(r,A,rho,C,D,gamma,r0)",
    "Potential-Form:density(r,n)",
    "EAM-Embed:Th",
    "EAM-Embed:U",
    "EAM-Embed:O",
    "EAM-Density:Th",
    "EAM-Density:U",
    "EAM-Density:O",
    "Pair:O-O",
    "Pair:Th-Th",
    "Pair:U-U",
    "Pair:Th-O",
    "Pair:U-O" ]

  expect.sort()
  
  with _get_lammps_resource_dir().join("CRG_U_Th.aspot").open() as infile:
    cp = ConfigParser(infile)
    items = _query_actions._list_item_labels(cp)
    items.sort()
    assert expect == items

def test_list_items():
  expect = [
    ("Tabulation:target", "setfl"),
    ("Tabulation:nr", "1000"),
    ("Tabulation:dr", "0.01"),
    ("Tabulation:nrho", "1000"),
    ("Tabulation:drho", "0.01"),
    ("Potential-Form:buck_morse(r,A,rho,C,D,gamma,r0)", "as.buck(r,A,rho,C) + as.morse(r, gamma, r0, D)" ),
    ("Potential-Form:density(r,n)", "(n/r^8) * 0.5 * (1+erf(20*(r-1.5)))" ),
    ("EAM-Embed:Th", "as.sqrt -1.185"),
    ("EAM-Embed:U", "as.sqrt -1.806"),
    ("EAM-Embed:O", "as.sqrt -0.690"),
    ("EAM-Density:Th", "density 1742.622"),
    ("EAM-Density:U", "density 3450.995"),
    ("EAM-Density:O", "density 106.856"),
    ("Pair:O-O","as.buck    830.283 0.352856 3.884372"),
    ("Pair:Th-Th","as.buck 18600 0.2884 0.0"),
    ("Pair:U-U","as.buck 18600 0.2747 0.0"),
    ("Pair:Th-O","buck_morse 315.544 0.395903 0.0 0.62614 1.85960 2.49788"),
    ("Pair:U-O" ,"buck_morse 448.779 0.387758 0.0 0.66080 2.05815 2.38051")]
  expect.sort()
  
  with _get_lammps_resource_dir().join("CRG_U_Th.aspot").open() as infile:
    cp = ConfigParser(infile)
    items = _query_actions._list_items(cp)
    items.sort()
    assert expect == items
    
def test_item_value():
  with _get_lammps_resource_dir().join("CRG_U_Th.aspot").open() as infile:
    cp = ConfigParser(infile)
    value = _query_actions._item_value(cp, "Tabulation:drho")
    assert "0.01" == value

def test_list_plot_item_labels():
  expect = [
    "EAM-Embed:Th",
    "EAM-Embed:U",
    "EAM-Embed:O",
    "EAM-Density:Th",
    "EAM-Density:U",
    "EAM-Density:O",
    "Pair:O-O",
    "Pair:Th-Th",
    "Pair:U-U",
    "Pair:Th-O",
    "Pair:U-O" ]

  expect.sort()
  
  with _get_lammps_resource_dir().join("CRG_U_Th.aspot").open() as infile:
    cp = ConfigParser(infile)
    items = _query_actions._list_plot_item_labels(cp)
    items.sort()
    assert expect == items

def test_key_normalisation():

  import io

  cfg1 = u"""[Potential-Form]
buck_morse(r, A,  rho,  C, D,  gamma,    r0) : test

[EAM-Density]
A -> B : test"""

  expect = [("Potential-Form:buck_morse(r,A,rho,C,D,gamma,r0)", "test"),
    ("EAM-Density:A->B", "test")
  ]

  expect.sort()

  with io.StringIO(cfg1) as infile:
    cp = ConfigParser(infile)
    items = _query_actions._list_items(cp)
    items.sort()
    assert expect == items