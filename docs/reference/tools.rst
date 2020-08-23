******************
Command Line Tools
******************

.. _potable-tool:

`potable`
=========

The :program:`potable` tool is the interface for working with potential definition files. In addition to converting a potential model definition into a tabulation it allows their contents to be queried, filtered, overridden and plotted.

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

.. option:: -h, --help

  show this help message and exit


.. rubric:: Query


Query items in the configuration file

.. option:: --list-items, -l

      List items in configuration file to STD_OUT. One is listed per line with format ``SECTION_NAME:KEY=VALUE``

.. option:: --list-item-labels

    List item in configuration file to STD_OUT. One item per line with format ``SECTION_NAME:KEY``


.. option:: --item-value SECTION_NAME:KEY

    Return the value for given item in configuration file

.. rubric:: Filter


Filter items from the configuration file


.. option:: --include-species [SPECIES [SPECIES ...]]

   If specified, only those ``SPECIES`` provided will be included in tabulation.

.. option:: --exclude-species [SPECIES [SPECIES ...]]

   ``SPECIES`` provided to this option will NOT be included in tabulation.


.. rubric:: Override


Add or override values in the configuration file

.. _cmdoption-override-item:

.. option:: --override-item [SECTION_NAME:KEY=VALUE [SECTION_NAME:KEY=VALUE ...]] , -e [SECTION_NAME:KEY=VALUE [SECTION_NAME:KEY=VALUE ...]]

  Use ``VALUE`` for item ``SECTION_NAME:KEY`` instead of value contained in the configuration file

.. option:: --add-item [SECTION_NAME:KEY=VALUE [SECTION_NAME:KEY=VALUE ...]] , -a [SECTION_NAME:KEY=VALUE [SECTION_NAME:KEY=VALUE ...]]

  Add item to configuration file 

.. option:: --remove-item [SECTION_NAME:KEY [SECTION_NAME:KEY ...]] , -r [SECTION_NAME:KEY [SECTION_NAME:KEY ...]]

   Remove item from configuration file


Examples:
---------

Various examples of the use of this tool are given throughout the documentation:

  * :ref:`quick-start`\ .
  * :ref:`Quick Start: Generating Basak Tabulation for DL_POLY <quickstart_basak_generating_tabulation>`\ .
  * :ref:`Quick Start: Generating Basak Tabulation for LAMMPS <specifying-other-tabulation-targets>`\ .
  * :ref:`Quick Start: Generating Basak Tabulation for GULP <specifying-other-tabulation-targets-GULP>`\ . This provides an example of the ``--override-item`` option.
  * :ref:`User Guide: Making and Testing the Tabulation - Sutton Ag Example <user_guide_potable_sutton_ag_making_and_testing>`\ .
  * :ref:`potable-troubleshooting`\ : provides an example of the ``--override-item`` option.
