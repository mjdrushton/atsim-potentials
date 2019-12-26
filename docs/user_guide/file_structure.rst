File Structure
--------------

.. toctree::

	pair_tabulation_section
	pair_pair_section
	pair_potential_form_section
	pair_table_form


Potable files for pair-potential models may contain the following sections::

	[Tabulation]
	...

	[Pair]
	...

	[Potential-Form]
	...

	[Table-Form:NAME]
	...

	[Species]
	...


Each section may contain a number of configuration options. These have the general form::

	ITEM : VALUE

or an equals sign may be used instead::

	ITEM = VALUE

Where ``ITEM`` identifies the configuration option's name and ``VALUE`` its value. 

Lines maybe commented out using ``#`` and line continuation is supported according to identation (see `here <https://docs.python.org/3/library/configparser.html#supported-ini-file-structure>`_ for more details). 


Each section of the file has a specific purpose and not all sections will be required in all cases:

* :ref:`potable_tabulation-section` 
	+ describes how the file should be converted into a table file by the :ref:`potable <potable-tool>` command. Contains information such as cutoff, output table format and cutoff. 
* :ref:`potable-pair-section` 
	+ this is where pair interactions are defined by parametrising a potential-form.
* :ref:`potable-potential-form` 
	+ this section allows custom potential-forms to be defined. This may be required when you can't find an appropriate function from those :ref:`supplied <list-of-potential-forms>` with ``atsim.potentials``\ . However in many cases this won't be necessary and this section needn't appear in your model definition.
* :ref:`potable-table-form` 
	+ Multiple ``[Table-Form]`` sections may be specified. These allow potential-forms to be defined from tables of x,y points. These tabulated forms can be used in the same way as potnetials defined in ``[Potential-Form]``\ .
* ``[Species]`` 
	+ this is used to provide meta-data about the species being tabulated. In most cases this section can be omitted (as very little species data is used during pair-tabulation). It is however, sometimes useful to include atomic charges etc, here so that the input file represents a complete description of a given potential model.

