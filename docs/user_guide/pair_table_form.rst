.. _potable-table-form:

************************
``[Table-Form]`` section
************************

This section of the input  serves a similar purpose to the :ref:`[Potential-Form] <potable-potential-form>` section as it allows custom potential functions to be defined. However, instead of being defined using a mathematical expression they are specified as tables of x,y points with interpolation providing intermediate values.

It is sometimes to convenient to use pre-tabulated values for very complex expressions however you should always use caution. It is recommended that when using ``[Table-Form]`` that you plot the resulting functions and derivatives. This is because interpolation can sometimes introduce spurious effects, so it's worth checking that nothing odd has happened.

.. _potable-potential-form-formula-syntax:

Defining a Table Form
+++++++++++++++++++++

The general form of a ``[Table-Form]`` section is::

	[Table-Form:NAME]
	interpolation : INTERPOLATION_TYPE
	xy : DATA

Multiple ``[Table-Form]`` sections can appear in a ``potable`` input file and are distinguished by the ``NAME`` identifier included after the colon in the section header. This identifier can be used elsewhere, such as in ``[Pair]`` to use the function.

The value of ``INTERPOLATION_TYPE`` specifies how values are calculated between data-points. A list of supported interpolation schemes can be found :ref:`here <ref-potable-input-table-form-interpolation>`\ , but for the examples that follow we will use ``cubic_spline`` interpolation.

Finally the function's data is provided through the ``xy`` option. The value of ``DATA`` is a list of space separated x,y pairs. To aid readability these can appear on separate lines as long as they are indented. Alternatively, table data can be specified as separate arrays of x and y values :ref:`see here for more <ref-potable-input-table-form-x>` \ . 




Example: [Table-Form] pair-potential
------------------------------------

Revisiting the example from earlier (see :ref:`basak-potential-model`\ ) the following shows how the O-U interaction from the Basak model [Basak2003]_ can be represented as a ``[Table-Form]``\ :


The data points from :ref:`fig_table_form_plot` (blue) have been included in the ``potable`` input file :download:`basak_table_form.aspot <example_files/basak_table_form.aspot>` using the :ref:`xy <ref-potable-input-table-form-xy>` option.

The third line of the ``[Pair]`` section deserves notice::

	O-U = tabulated

The potential-form ``tabulated`` specified for the ``O-U`` interaction refers to the name of the table-form: ``[Table-Form:tabulated]``\ . It should also be noted that instances of the potential form do not take any parameters.

.. literalinclude:: example_files/basak_table_form.aspot


The following figure shows the good match between the analytical and tabulated forms resulting from the use of the ``[Table-Form]``\ .

.. figure:: figures/table_form_plot.*
	:name: fig_table_form_plot

	Plot showing the table points used in the example (blue circles) and the function resulting from cubic spline interpolation (red). The original analytical form is shown in grey.