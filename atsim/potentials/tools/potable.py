"""Front-end script for atsim.potentials. 
Allows potentials to be tabulated using simple .ini based configuration files"""

import argparse
import itertools

import logging

from ..config import Configuration, ConfigParser

def _list_section(cp, section):
  outlist = []
  raw_cp = cp.raw_config_parser
  for k in raw_cp[section]:
    v = raw_cp[section][k]
    k = "{section}:{key}".format(section = section, key = k)
    ov = (k,v)
    outlist.append(ov)
  return outlist

def _list_pair(cp):
  return _list_section(cp, "Pair")

def _list_potential_form(cp):
  return _list_section(cp, "Potential-Form")

def _list_tabulation(cp):
  return _list_section(cp, "Tabulation")

def _list_eam_dens(cp):
  return _list_section(cp, "EAM-Density")

def _list_eam_embed(cp):
  return _list_section(cp, "EAM-Embed")

def _parse_raw(cp, orphan_sections):
  outlist = []
  for s in orphan_sections:
    entries = _list_section(cp, s)
    outlist.extend(entries)
  return outlist

def _list_items(cfg_file):
  cp = ConfigParser(cfg_file)
  items = []
  parsed_sections = cp.parsed_sections

  if "pair" in parsed_sections:
    pair_items = _list_pair(cp)
    items.extend(pair_items)
  
  if "potential_form" in parsed_sections:
    pf_items = _list_potential_form(cp)
    items.extend(pf_items)

  if "tabulation" in parsed_sections:
    tab_items = _list_tabulation(cp)
    items.extend(tab_items)
  
  if "eam_embed" in parsed_sections:
    embed_items = _list_eam_embed(cp)
    items.extend(embed_items)

  if ("eam_density" in parsed_sections) or ("eam_density_fs" in parsed_sections):
    eam_dens_items = _list_eam_dens(cp)
    items.extend(eam_dens_items)

  orphan_sections = cp.orphan_sections
  raw_items = _parse_raw(cp, orphan_sections)
  items.extend(raw_items)
  return items

def _list_item_labels(cfg_file):
  items = _list_items(cfg_file)
  outlist = [k for (k,v) in items]
  return outlist

def _list_plot_item_labels(cfg_file):
  cp = ConfigParser(cfg_file)
  items = list(itertools.chain(
    _list_pair(cp),
    _list_eam_dens(cp),
    _list_eam_embed(cp)
  ))
  outlist = [k for (k,v) in items]
  return outlist  

def _item_value(cp, key):
  section, section_key = key.split(":",1)
  v = cp.raw_config_parser[section][section_key]
  return v 

def _parse_command_line():
  p = argparse.ArgumentParser(description = "Tabulate potential models for common atomistic simulation codes. This is part of the atsim.potentials package.")

  p.add_argument("in_filename", metavar = "POTENTIAL_DEFN_FILE", help = "File containing definition of potential model.")
  p.add_argument("out_filename", metavar = "OUTPUT_FILE", help = "File into which data will be tabulated.")
  # TODO: Provide options to override the contents of the [Tabulation] section of config file.

  args = p.parse_args()
  return args

def _setup_logging():
  logging.basicConfig(level = logging.INFO, format = "%(message)s")

def main():
  _setup_logging()
  logger = logging.getLogger(__name__).getChild("main")
  args = _parse_command_line()

  ini_filename = args.in_filename
  out_filename = args.out_filename

  config = Configuration()
  with open(ini_filename) as infile:
    tabulation = config.read(infile)
  
  with open(out_filename, 'w') as outfile:
    logger.info("Writing output to: {}".format(out_filename))
    tabulation.write(outfile)

if __name__ == '__main__':
  main()