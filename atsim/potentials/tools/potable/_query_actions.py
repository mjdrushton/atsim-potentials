import itertools
from ...config import Configuration, ConfigParser

import sys

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

def _list_items(cp):
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

def _list_item_labels(cp):
  items = _list_items(cp)
  outlist = [k for (k,v) in items]
  return outlist

def _list_plot_item_labels(cp):
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

def action_list_items(cp):
  items = _list_items(cp)
  for k, v in items:
    outline = "{}={}\n".format(k,v)
    sys.stdout.write(outline)

def action_list_item_labels(cp):
  items = _list_items(cp)
  for k, v in items:
    outline = "{}\n".format(k)
    sys.stdout.write(outline)
  
def action_item_value(cp, key):
  v = _item_value(cp, key)
  outline = "{}\n".format(v)
  sys.stdout.write(outline)