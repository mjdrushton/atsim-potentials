
from setuptools import setup, find_packages

setup(name="atsim.potentials",
  use_hg_version=True,

  packages = find_packages(exclude=["tests"]),
  namespace_packages = ["atsim"],
  test_suite = "tests",
  setup_requires = ["hgtools"],
  install_requires = ["setuptools"],

  #zip_safe = True,

  # Meta-data for PyPI
  description = "atsim.potentials provides tools for working with pair and embedded atom method potential models including tabulation routines for DL_POLY and LAMMPS",
  long_description = open('README.rst').read(),
  author = "M.J.D. Rushton",
  author_email = "m.rushton@imperial.ac.uk",
  license = "Apache License (2.0)",
  url = "https://bitbucket.org/mjdr/atsim_potentials",
  download_url = "https://bitbucket.org/mjdr/atsim_potentials/get/0.0.0.tar.gz",
  keywords = [
  "pair potentials",
    "embedded atom model",
    "LAMMPS",
    "DL_POLY",
    "potential tabulation",
    "atomistic",
    "simulation",
    "molecular dynamics",
    "atomic scale"],
  classifiers = [
    "License :: OSI Approved :: Apache Software License",
    "Programming Language :: Python :: 2.7",
    "Topic :: Scientific/Engineering"])
