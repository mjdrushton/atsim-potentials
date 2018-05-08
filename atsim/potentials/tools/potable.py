"""Front-end script for atsim.potentials. 
Allows potentials to be tabulated using simple .ini based configuration files"""

import argparse

import logging

from ..config import Configuration

def _parse_command_line():
  p = argparse.ArgumentParser(description = "Tabulate potential models for common atomistic simulation codes. This is part of the atsim.potentials package.")

  p.add_argument("in_filename", metavar = "POTENTIAL_DEFN_FILE", help = "File containing definition of potential model.")
  p.add_argument("out_filename", metavar = "OUTPUT_FILE", help = "File into which data will be tabulated.")
  # TODO: Provide options to override the contents of the [Tabulation] section of config file.

  args = p.parse_args()

  return args


def _setup_logging():
  logging.basicConfig(level = logging.INFO)

def main():
  _setup_logging()
  args = _parse_command_line()

  ini_filename = args.in_filename
  out_filename = args.out_filename

  config = Configuration()
  with open(ini_filename) as infile:
    tabulation = config.read(infile)
  
  with open(out_filename, 'w') as outfile:
    tabulation.write(outfile)

if __name__ == '__main__':
  main()