
from setuptools import setup, find_packages

setup(name="atsim.potentials",
  packages = find_packages(exclude=["tests"]),
  namespace_packages = ["atsim"],
  test_suite = "tests",
  install_requires = [
    "setuptools",
    "configparser",
    "funcsigs",
    "pyparsing>=2.2.0",
    "cexprtk>=0.3.4",
    "wrapt==1.11.2",
    "scipy"],

  version = '0.3.0dev',

  entry_points = {
        'console_scripts' : [
          'potable=atsim.potentials.tools.potable:main'
        ]
  },

  # Meta-data for PyPI
  description = "atsim.potentials provides tools for working with pair and embedded atom method potential models including tabulation routines for DL_POLY and LAMMPS",
  long_description = open('README.rst').read(),
  author = "M.J.D. Rushton",
  author_email = "m.j.d.rushton@gmail.com",
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
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.4",
    "Programming Language :: Python :: 3.6",
    "Programming Language :: Python :: 3.7",
    "Programming Language :: Python :: 3.8",
    "Topic :: Scientific/Engineering"])
