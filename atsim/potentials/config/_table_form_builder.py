import inspect

from  .. import tableforms
from ..potentialforms import is_potential
from ._common import Table_Form_Exception
from ._potential_form import Existing_Potential_Form
from ._python_potential_function import _Python_Potential_Function

class Table_Form(Existing_Potential_Form):

  def __init__(self, table_form_tuple, table_form_cls):
    self._pf = Table_Form_Factory(table_form_tuple, table_form_cls)
    Existing_Potential_Form.__init__(self, table_form_tuple.name, self._pf)
    self._func = _Python_Potential_Function(self.potential_definition, self._pf.potential_function)

  # @property
  # def register_function(self, func):
  #   pass

  @property
  def potential_function(self):
    return self._func

class Table_Form_Factory(object):
  """Instantiates and holds a Table_Form instance which is returned when this object's call
  method is invoked"""

  def __init__(self, table_form_tuple, table_form_cls):
    self._obj = table_form_cls(table_form_tuple.x, table_form_tuple.y)

  def __call__(self):
    return self._obj

  @property
  def potential_function(self):
    return self._obj
    

class Table_Form_Builder(object):
  """Class used by Potential_Form_Registry to create Potential_Form objects
  from the TableFormTuple instances returned by ConfigParser.table_form """

  def __init__(self):
    # Create a mapping from interpolation labels to class
    self._table_forms = self._populate()

  def _populate(self):
    tf_dict = {}
    classes  = inspect.getmembers(tableforms, is_potential)
    for cls_name, cls in classes:
      label = cls.config_label
      tf_dict[label] = cls
    return tf_dict

  def _config_name_to_class(self, interpolation_type):
    """Return table form class for a given interpolation_type"""
    return self._table_forms[interpolation_type]

  def create_potential_form(self, table_tuple):
    """Return an instance of Existing_Potential_Form based on table form definition
    given by `table_tuple`.

    :param table_tuple: TableFormTuple instance defining table form.
    :return: Potential_Form instance"""
    name = table_tuple.name
    try:
      cls = self._config_name_to_class(table_tuple.interpolation)
    except KeyError:
      raise Table_Form_Exception("Unknown interpolation type specified for [Table-Form:{}]: '{}'".format(table_tuple.name, table_tuple.interpolation))

    factory = Table_Form_Factory(table_tuple, cls)
    func = factory()
    pf = Table_Form(table_tuple, cls)
    return pf