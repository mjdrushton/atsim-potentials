
from setuptools import setup, find_packages

setup(name="atsim.potentials",
  packages = find_packages(exclude=["tests"]),
  namespace_packages = ["atsim"],
  test_suite = "tests",
  install_requires = [
    'funcsigs==1.0.2;python_version<"3.5"',
    'openpyxl==2.6.4;python_version=="3.5.*"',
    'openpyxl>=3.0.0;python_version>="3.6"',
    "cexprtk>=0.3.4",
    "pyparsing>=2.2.0",
    "scipy",
    "setuptools",
    "wrapt>=1.12.1",
    ],

  version = '0.4.1',

  entry_points = {
        'console_scripts' : [
          'potable=atsim.potentials.tools.potable:main'
        ]
  },

  # Meta-data for PyPI
  description = "atsim.potentials provides tools for working with pair and embedded atom method potential models including tabulation routines for DL_POLY and LAMMPS",
  long_description_content_type= 'text/x-rst',
  long_description = open('README.rst').read(),
  author = "M.J.D. Rushton",
  author_email = "m.j.d.rushton@gmail.com",
  license = "Apache License (2.0)",
  url = "https://github.com/mjdrushton/atsim-potentials",
  download_url = "https://github.com/mjdrushton/atsim-potentials/archive/master.zip",
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
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.6",
    "Programming Language :: Python :: 3.7",
    "Programming Language :: Python :: 3.8",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3.10",
    "Topic :: Scientific/Engineering"])
