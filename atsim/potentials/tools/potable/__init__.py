"""Front-end script for atsim.potentials. 
Allows potentials to be tabulated using simple .ini based configuration files"""

import argparse
import logging

import sys

from . import _query_actions
from . import _actions

def _parse_command_line():
  p = argparse.ArgumentParser(description = "Tabulate potential models for common atomistic simulation codes. This is part of the atsim.potentials package.")

  p.add_argument("config_file", metavar = "POTENTIAL_DEFN_FILE", type=argparse.FileType('r'), help = "File containing definition of potential model.")
  p.add_argument("out_filename", nargs = '?', metavar = "OUTPUT_FILE", help = "File into which data will be tabulated.")
  # TODO: Provide options to override the contents of the [Tabulation] section of config file.

  query_group = p.add_argument_group("Query", "Query items in the configuration file")
  qmutex_group = query_group.add_mutually_exclusive_group()
  qmutex_group.add_argument("--list-items", action = 'store_true', help = "List items in configuration file to STD_OUT. One is listed per line with format SECTION_NAME:KEY=VALUE")
  qmutex_group.add_argument("--list-item-labels", action = 'store_true', help = "List item in configuration file to STD_OUT. One item per line with format SECTION_NAME:KEY")
  qmutex_group.add_argument("--item-value",  type = str, nargs = 1, metavar = "SECTION_NAME:KEY", help = "Return the value for given item in configuration file")

  args = p.parse_args()
  return p, args

def _setup_logging():
  logging.basicConfig(level = logging.INFO, format = "%(message)s")

def main():
  _setup_logging()
  logger = logging.getLogger(__name__).getChild("main")
  p, args = _parse_command_line()

  config_file = args.config_file
  
  if args.list_items:
    _query_actions.action_list_items(config_file)
    sys.exit(0)
  elif args.list_item_labels:
    _query_actions.action_list_item_labels(config_file)
    sys.exit(0)
  elif args.item_value:
    _query_actions.action_item_value(config_file, args.item_value[0])
    sys.exit(0)


  # Otherwise just perform tabulation
  if not args.out_filename:
    p.error("Path of OUTPUT_FILE for tabulation not specified.")
  _actions.action_tabulate(config_file, args.out_filename)
  sys.exit(0)

if __name__ == '__main__':
  main()