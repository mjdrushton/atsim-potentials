# Change Log

## 0.4.1 (2022-12-22)

### Bug-Fixes

* Fixes to ``setup.py`` to allow atsim.potentials to install under python 3.10
* Test for factorial function in maths expressions was failing under python 3.10. This is now fixed.

## 0.4.0 (2021-5-20)

### New Features

* Added support for angular dependent potential (ADP) models (LAMMPS `pair_style adp`).
* Support for `[Variables]` section in `potable` files which allows use of string snippets and string interpolation.

## 0.3.0 (2020-8-24)
### New Features

* Introduced the new potable tool to allow tabulation without needing to write a python script.
* Revamped python api.

### Bug-Fixs

* The order that the density functions were specified to the writeSetFLFinnisSinclair() function was the reverse of what would be expected. This has been fixed.

## 0.2.1 (2018-05-19)
### Bug-Fixes

* Fix to the plot functions. Previously all x-axis values were being set to the `lowx` value given to the function.

## 0.2.0 (2018-05-01)
### New Features

* Support for Python 3

## 0.1.1 (2014-03-25)

* Initial Release
