******************
Command Line Tools
******************

.. _potable-tool:

`potable`
=========

The `potable` tool is the interface for working with potential definition files. In addition to converting a potential model definition into a tabulation it allows their contents to be queried, filtered, overridden and plotted.

Usage
-----

::

	potable [-h]
	               [--list-items | --list-item-labels | --item-value SECTION_NAME:KEY]
	               [--include-species [SPECIES [SPECIES ...]] | --exclude-species
	               [SPECIES [SPECIES ...]]]
	               [--override-item [SECTION_NAME:KEY=VALUE [SECTION_NAME:KEY=VALUE ...]]]
	               [--add-item [SECTION_NAME:KEY=VALUE [SECTION_NAME:KEY=VALUE ...]]]
	               [--remove-item [SECTION_NAME:KEY [SECTION_NAME:KEY ...]]]
	               POTENTIAL_DEFN_FILE [OUTPUT_FILE]


Tabulate potential models for common atomistic simulation codes. This is part
of the atsim.potentials package.

Positional Arguments:
~~~~~~~~~~~~~~~~~~~~~

  * ``POTENTIAL_DEFN_FILE``   File containing definition of potential model.
  * ``OUTPUT_FILE``           File into which data will be tabulated.


Optional Arguments:
~~~~~~~~~~~~~~~~~~~

  * ``-h, --help``            show this help message and exit

Query:
++++++

Query items in the configuration file

  * ``--list-items``, ``-l``      List items in configuration file to STD_OUT. One is
                        listed per line with format ``SECTION_NAME:KEY=VALUE``
  * ``--list-item-labels``    List item in configuration file to STD_OUT. One item
                        per line with format ``SECTION_NAME:KEY``
  * ``--item-value SECTION_NAME:KEY``
                        Return the value for given item in configuration file

Filter:
+++++++

Filter items from the configuration file

  * ``--include-species [SPECIES [SPECIES ...]]`` If specified, only those ``SPECIES`` provided will be included in tabulation.
  * ``--exclude-species [SPECIES [SPECIES ...]]`` ``SPECIES`` provided to this option will NOT be included in tabulation.

Override:
+++++++++

Add or override values in the configuration file

  * ``--override-item [SECTION_NAME:KEY=VALUE [SECTION_NAME:KEY=VALUE ...]]``\ , ``-e [SECTION_NAME:KEY=VALUE [SECTION_NAME:KEY=VALUE ...]]`` Use ``VALUE`` for item ``SECTION_NAME:KEY`` instead of value contained in the configuration file
  * ``--add-item [SECTION_NAME:KEY=VALUE [SECTION_NAME:KEY=VALUE ...]]``\ , ``-a [SECTION_NAME:KEY=VALUE [SECTION_NAME:KEY=VALUE ...]]`` Add item to configuration file 
  * ``--remove-item [SECTION_NAME:KEY [SECTION_NAME:KEY ...]]``\ , ``-r [SECTION_NAME:KEY [SECTION_NAME:KEY ...]]`` Remove item from configuration file


Examples:
---------

.. todo::
	
	Provide examples of potable's usage.