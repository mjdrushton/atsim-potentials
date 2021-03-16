import cexprtk

from ._common import Potential_Form_Exception


class _Cexptrk_Potential_Function(object):
  """Callable that can be added to cexprtk symbol_table. 

  Its wrapped cexprtk expression is only instantiated on first use. This is to allow
  the expression's symbol table to be populated before the expression is parsed. 
  Otherwise custom functions may not have been defined at the time of parsing."""

  def __init__(self, potential_form_tuple):
    """:param potential_form_tuple: PotentialFormTuple describing this function"""
    self._potential_form_tuple = potential_form_tuple
    self._local_symbol_table = self._init_symbol_table()
    self._expression = None

  def _init_symbol_table(self):
    local_symbol_table = cexprtk.Symbol_Table({}, add_constants = True)
    parameter_names = self._potential_form_tuple.signature.parameter_names
    for pn in parameter_names:
      local_symbol_table.variables[pn] = 1.0
    return local_symbol_table

  def register_function(self, func):
    """Register `func` with this object's symbol_table"""
    label = func._potential_form_tuple.signature.label
    try:
      self._local_symbol_table.functions[label] = func
    except cexprtk._exceptions.NameShadowException as e:
      msg = "Name clash for potential-form '{}': {}".format(label, str(e))
      raise Potential_Form_Exception(msg)
      

  def __call__(self, *args):
    parameter_names = self._potential_form_tuple.signature.parameter_names
    assert len(args) == len(parameter_names)
    for (pn, v) in zip(parameter_names, args):
      self._local_symbol_table.variables[pn] = v

    try:
      if not self._expression:
        try:
          self._expression = cexprtk.Expression(self._potential_form_tuple.expression, self._local_symbol_table)
        except cexprtk.ParseException as pe:
          raise Potential_Form_Exception("mathematical expression couldn't be parsed {}".format(pe))
      retval = self._expression()
      return retval
    except Potential_Form_Exception as e:
      msg = e.args[0]
      sig = ",".join(self._potential_form_tuple.signature.parameter_names)
      sig = "{label}({sig})".format(label = self._potential_form_tuple.signature.label, sig = sig)
      msg = "In potential-form '{sig} = {expression}': {msg}".format(
        msg = msg, 
        sig = sig, 
        expression = self._potential_form_tuple.expression)
      raise Potential_Form_Exception(msg)
      



