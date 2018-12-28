"""Front-end script for atsim.potentials. 
Allows potentials to be tabulated using simple .ini based configuration files"""

import argparse
import logging
import sys
import collections
import itertools

from . import _query_actions
from . import _actions

from ...config import ConfigParserOverrideTuple, ConfigParser, FilteredConfigParser
from ...config._common import ConfigurationException

def _parse_command_line(cli_args = None):
  p = argparse.ArgumentParser(description = "Tabulate potential models for common atomistic simulation codes. This is part of the atsim.potentials package.")

  p.add_argument("config_file", metavar = "POTENTIAL_DEFN_FILE", type=argparse.FileType('r'), help = "File containing definition of potential model.")
  p.add_argument("out_filename", nargs = '?', metavar = "OUTPUT_FILE", help = "File into which data will be tabulated.")

  query_group = p.add_argument_group("Query", "Query items in the configuration file")
  qmutex_group = query_group.add_mutually_exclusive_group()
  qmutex_group.add_argument("--list-items", "-l", action = 'store_true', help = "List items in configuration file to STD_OUT. One is listed per line with format SECTION_NAME:KEY=VALUE")
  qmutex_group.add_argument("--list-item-labels", action = 'store_true', help = "List item in configuration file to STD_OUT. One item per line with format SECTION_NAME:KEY")
  qmutex_group.add_argument("--item-value",  type = str, nargs = 1, metavar = "SECTION_NAME:KEY", help = "Return the value for given item in configuration file")
  
  filter_group = p.add_argument_group("Filter", "Filter items from the configuration file")
  fmutex_group = filter_group.add_mutually_exclusive_group()
  fmutex_group.add_argument("--include-species", nargs = '*', metavar = "SPECIES", help = "If specified, only those SPECIES provided will be included in tabulation.")
  fmutex_group.add_argument("--exclude-species", nargs = '*', metavar = "SPECIES", help = "SPECIES given provided to this option will NOT be included in tabulation.")

  override_group = p.add_argument_group("Override", "Add or override values in the configuration file")
  override_group.add_argument("--override-item", "-e", nargs='*', action="append", metavar = "SECTION_NAME:KEY=VALUE", help = "Use VALUE for item SECTION_NAME:KEY instead of value contained in the configuration file")
  override_group.add_argument("--add-item", "-a", nargs='*', action="append", metavar = "SECTION_NAME:KEY=VALUE", help = "Add item to configuration file")
  override_group.add_argument("--remove-item", "-r", nargs='*', action="append", metavar = "SECTION_NAME:KEY", help = "Remove item from configuration file")

  args = p.parse_args(args = cli_args)
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

def _make_config_parser(cfg_file, overrides, additional, remove, species, exclude_flag):
  override_dict = collections.OrderedDict()
  if not overrides is None:
    for override in itertools.chain.from_iterable(overrides):
      over_tuple = _create_override_tuple(override)
      k = (over_tuple.section, over_tuple.key)
      override_dict[k] = over_tuple

  if not remove is None:
    for override in itertools.chain.from_iterable(remove):
      over_tuple = _create_override_tuple(override, False)
      k = (over_tuple.section, over_tuple.key)
      override_dict[k] = over_tuple

  overrides_list = list(override_dict.values())

  additional_list = []
  if not additional is None:
    for add in itertools.chain.from_iterable(additional):
      over_tuple = _create_override_tuple(add)
      additional_list.append(over_tuple)

  cp = ConfigParser(cfg_file, overrides = overrides_list, additional= additional_list)

  if species:
    if exclude_flag:
      cp = FilteredConfigParser(cp, exclude = species)
    else:
      cp = FilteredConfigParser(cp, include = species)
  return cp

def _setup_logging():
  logging.basicConfig(level = logging.INFO, format = "%(message)s")

def _do_tabulation(p, args):
  logger = logging.getLogger(__name__).getChild("main")
  species_list = None
  exclude_flag = False
  if args.include_species:
    species_list = args.include_species
  elif args.exclude_species:
    species_list = args.exclude_species
    exclude_flag = True

  cp = _make_config_parser(
    args.config_file, 
    args.override_item, 
    args.add_item, 
    args.remove_item,
    species_list, exclude_flag)
  
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

def main():
  _setup_logging()
  logger = logging.getLogger(__name__).getChild("main")
  p, args = _parse_command_line()

  try:
    _do_tabulation(p, args)
  except ConfigurationException as e:
    p.error("configuration error - {}".format(e))


if __name__ == '__main__':
  main()