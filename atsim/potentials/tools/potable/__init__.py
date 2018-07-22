"""Front-end script for atsim.potentials. 
Allows potentials to be tabulated using simple .ini based configuration files"""

import argparse
import logging
import sys
import collections

from . import _query_actions
from . import _actions

from ...config import ConfigParserOverrideTuple, ConfigParser

def _parse_command_line():
  p = argparse.ArgumentParser(description = "Tabulate potential models for common atomistic simulation codes. This is part of the atsim.potentials package.")

  p.add_argument("config_file", metavar = "POTENTIAL_DEFN_FILE", type=argparse.FileType('r'), help = "File containing definition of potential model.")
  p.add_argument("out_filename", nargs = '?', metavar = "OUTPUT_FILE", help = "File into which data will be tabulated.")

  query_group = p.add_argument_group("Query", "Query items in the configuration file")
  qmutex_group = query_group.add_mutually_exclusive_group()
  qmutex_group.add_argument("--list-items", "-l", action = 'store_true', help = "List items in configuration file to STD_OUT. One is listed per line with format SECTION_NAME:KEY=VALUE")
  qmutex_group.add_argument("--list-item-labels", action = 'store_true', help = "List item in configuration file to STD_OUT. One item per line with format SECTION_NAME:KEY")
  qmutex_group.add_argument("--item-value",  type = str, nargs = 1, metavar = "SECTION_NAME:KEY", help = "Return the value for given item in configuration file")
  
  override_group = p.add_argument_group("Override", "Add or override values in the configuration file")
  override_group.add_argument("--override-item", "-e", nargs='*', metavar = "SECTION_NAME:KEY=VALUE", help = "Use VALUE for item SECTION_NAME:KEY instead of the in the configuration file")
  override_group.add_argument("--add-item", "-a", nargs='*', metavar = "SECTION_NAME:KEY=VALUE", help = "Add item to configuration file")
  override_group.add_argument("--remove-item", "-r", nargs='*', metavar = "SECTION_NAME:KEY", help = "Remove item from configuration file")

  args = p.parse_args()
  return p, args

def _create_override_tuple(key, has_value = True):
  # TODO: Error handling for malformed options
  section,key = key.split(":", 1)
  if has_value:
    key, value = key.split("=", 1)
  else:
    value = None
  retval = ConfigParserOverrideTuple(section = section, key = key, value = value)
  return retval

def _make_config_parser(cfg_file, overrides, additional, remove):
  override_dict = collections.OrderedDict()
  if not overrides is None:
    for override in overrides:
      over_tuple = _create_override_tuple(override)
      k = (over_tuple.section, over_tuple.key)
      override_dict[k] = over_tuple

  if not remove is None:
    for override in remove:
      over_tuple = _create_override_tuple(override, False)
      k = (over_tuple.section, over_tuple.key)
      override_dict[k] = over_tuple

  overrides_list = list(override_dict.values())

  additional_list = []
  if not additional is None:
    for add in additional:
      over_tuple = _create_override_tuple(add)
      additional_list.append(over_tuple)

  cp = ConfigParser(cfg_file, overrides = overrides_list, additional= additional_list)
  return cp

def _setup_logging():
  logging.basicConfig(level = logging.INFO, format = "%(message)s")

def main():
  _setup_logging()
  logger = logging.getLogger(__name__).getChild("main")
  p, args = _parse_command_line()
  cp = _make_config_parser(args.config_file, args.override_item, args.add_item, args.remove_item)
  
  if args.list_items:
    _query_actions.action_list_items(cp)
    sys.exit(0)
  elif args.list_item_labels:
    _query_actions.action_list_item_labels(cp)
    sys.exit(0)
  elif args.item_value:
    _query_actions.action_item_value(cp, args.item_value[0])
    sys.exit(0)

  # Otherwise just perform tabulation
  if not args.out_filename:
    p.error("Path of OUTPUT_FILE for tabulation not specified.")
  _actions.action_tabulate(cp, args.out_filename)
  sys.exit(0)

if __name__ == '__main__':
  main()