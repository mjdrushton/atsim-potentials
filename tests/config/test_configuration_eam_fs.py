"""Tests for Finnis-Sinclair style potentials"""

from atsim.potentials.config import Configuration
from atsim.potentials.config._config_parser import ConfigParser
from atsim.potentials.config._configuration import Configuration
from atsim.potentials.potentialforms import zero

from .test_configuration_eam import _get_lammps_resource_dir

from .._eam_fs_AlFe import alEmbedFunction
from .._eam_fs_AlFe import feEmbedFunction

from .._eam_fs_AlFe import alAlDensFunction
from .._eam_fs_AlFe import feFeDensFunction
from .._eam_fs_AlFe import feAlDensFunction

from .._eam_fs_AlFe import ppfuncAlAl, ppfuncAlFe, ppfuncFeFe

import pytest
import py

import itertools

def _get_resource_dir():
  rd = py.path.local(__file__).dirpath()
  rd = rd.join('config_resources')
  return rd


def test_fs_AlFe():
  cfg_file_path = _get_lammps_resource_dir().join("AlFe_setfl_fs.aspot")
  with cfg_file_path.open() as config_file:
    config_parser = ConfigParser(config_file)
  assert config_parser.tabulation.target == "setfl_fs"

  with cfg_file_path.open() as config_file:
    cfgobj = Configuration()
    tabulation = cfgobj.read(config_file)

  # Now check that the EAM potentials are as they should be.
  assert tabulation.target == "setfl_fs"
  assert 2 == len(tabulation.eam_potentials)

  epots = tabulation.eam_potentials
  epot_dict = dict(zip([p.species for p in epots], epots))

  assert ["Al", "Fe"] == sorted([p.species for p in epots])  

  # Check the embedding functions
  alEmbedFunction_actual = epot_dict["Al"].embeddingFunction
  for rho in range(1, 300):
    expect = alEmbedFunction(rho)
    actual = alEmbedFunction_actual(rho)
    assert pytest.approx(expect) == actual

  feEmbedFunction_actual = epot_dict["Fe"].embeddingFunction
  for rho in range(1, 300):
    expect = feEmbedFunction(rho)
    actual = feEmbedFunction_actual(rho)
    assert pytest.approx(expect) == actual
    # print("{},{}".format(expect, actual))
  # pytest.fail()

  # Check the density functions
  al_pot  = epot_dict["Al"]
  dens_dict = al_pot.electronDensityFunction
  assert ["Al", "Fe"] == sorted(dens_dict.keys())

  alAlDensFunction_actual = dens_dict["Al"]
  alFeDensFunction_actual = dens_dict["Fe"]
  for r in range(65):
    r = r * 0.1
    expect = alAlDensFunction(r)
    actual = alAlDensFunction_actual(r)
    assert pytest.approx(expect) == actual

    expect = feAlDensFunction(r)
    actual = alFeDensFunction_actual(r)
    assert pytest.approx(expect) == actual

  fe_pot  = epot_dict["Fe"]
  dens_dict = fe_pot.electronDensityFunction
  assert ["Al", "Fe"] == sorted(dens_dict.keys())

  feAlDensFunction_actual = dens_dict["Al"]
  feFeDensFunction_actual = dens_dict["Fe"]
  for r in range(65):
    r = r * 0.1
    expect = feAlDensFunction(r)
    actual = feAlDensFunction_actual(r)
    assert pytest.approx(expect) == actual

    expect = feFeDensFunction(r)
    actual = feFeDensFunction_actual(r)
    assert pytest.approx(expect) == actual

  # Check pair-potentials
  ppots = tabulation.potentials

  pairs = [(ppot.speciesA, ppot.speciesB) for ppot in ppots]
  assert [("Al", "Al"), ("Al", "Fe"), ("Fe", "Fe")] == pairs

  expect_dict = { 
    ("Al", "Al") : ppfuncAlAl, 
    ("Al", "Fe") : ppfuncAlFe,
    ("Fe", "Fe") : ppfuncFeFe
  }

  for actual_pp in ppots:
    expect_pp = expect_dict[(actual_pp.speciesA, actual_pp.speciesB)]
    for r in range(100):
      r = r / 10.0
      expect = expect_pp(r)
      actual = actual_pp.energy(r)
      assert pytest.approx(expect) == actual

def test_add_zeroed_interactions_fs():
  with _get_resource_dir().join("spinel.aspot").open() as infile:
    cfg = Configuration()
    tab = cfg.read(infile)

    # Check that eam potential exist
    # expect_embed = set(['Ga', 'In', 'O', 'Mg', 'Al'])
    expect_embed = set(['Ga', 'In', 'O'])
    actual_embed = set([ep.species for ep in tab.eam_potentials])
    assert expect_embed == actual_embed

    # Check that the defaulted embedding functions are instances of as.zero
    zero_type = type(zero().func)
    as_zero_func_species = set([ep.species for ep in tab.eam_potentials if hasattr(ep.embeddingFunction, "func") and type(ep.embeddingFunction.func) == zero_type ])
    # assert set(['O', 'Mg', 'Al']) == as_zero_func_species
    assert set(['O']) == as_zero_func_species

    # Check that null-density functions exist.
    expect_density = set(itertools.product(expect_embed, expect_embed))

    actual_density = set()
    for ep in tab.eam_potentials:
      a = ep.species
      for b, func in ep.electronDensityFunction.items():
        actual_density.add( (a,b))
    
    assert expect_density == actual_density

    # Check that tabulation can complete.
    import io
    outfile = io.StringIO()
    tab.write(outfile)

    
def test_add_zeroed_interactions():
  with _get_resource_dir().join("setfl.aspot").open() as infile:
    cfg = Configuration()
    tab = cfg.read(infile)

    # Check that eam potential exist
    expect_embed = set(['Ga', 'In', 'O'])
    actual_embed = set([ep.species for ep in tab.eam_potentials])
    assert expect_embed == actual_embed

    # Check that the defaulted embedding functions are instances of as.zero
    zero_type = type(zero().func)
    as_zero_func_species = set([ep.species for ep in tab.eam_potentials if hasattr(ep.embeddingFunction, "func") and type(ep.embeddingFunction.func) == zero_type ])
    # assert set(['O', 'Mg', 'Al']) == as_zero_func_species
    assert set(['O']) == as_zero_func_species

    # Check that null-density functions exist.
    expect_density = expect_embed

    actual_density = set()
    for ep in tab.eam_potentials:
      a = ep.species
      actual_density.add(a)
    assert expect_density == actual_density

    # Check that tabulation can complete.
    import io
    outfile = io.StringIO()
    tab.write(outfile)
