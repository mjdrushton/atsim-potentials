"""Module defining the pyparsing grammar for multi-range and modified
potential definitions used by the config_parser."""

from pyparsing import *

def _grammar():
  from pyparsing import pyparsing_common
  number = pyparsing_common.number
  identifier = Combine(pyparsing_common.identifier+ZeroOrMore(Literal(".")+pyparsing_common.identifier))

  # multi_range
  range_start = Group((Literal(u">=") | Literal(u">"))("range_type") + number("start"))("range_start")
  # range_start = Group((unicodeString(u">=") | unicodeString(u">"))("range_type") + number("start"))("range_start")

  modified = Forward()
  potential_description = Group(identifier("potential_label") + Group(ZeroOrMore(number))("potential_parameters"))("potential_description")
  potential_definition =  modified | potential_description

  multi_range = Group(Optional(range_start) + potential_definition + ZeroOrMore(range_start + potential_definition))("multi_range")

  lpar = Literal('(').suppress()
  rpar = Literal(')').suppress()

  modifier_parameters = delimitedList(multi_range)

  modified << Group(identifier("modifier_label") + lpar +  Group(modifier_parameters)("modifier_parameters") + rpar)("modifier")
  # multi_range.setDebug()
  return  multi_range

multi_range_parser = _grammar()