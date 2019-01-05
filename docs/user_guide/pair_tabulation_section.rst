.. _potable_tabulation-section:

************************
``[Tabulation]`` section
************************

This section defines the file format in which ``potable`` will write its output files through the :ref:`target <ref-potable-input-tabulation-target>` configuration item.

.. seealso::

	* Input format reference: :ref:`ref-potable-input-tabulation`\ .



The pair-tabulation :ref:`target <ref-potable-input-tabulation-target>` can be one of:

* ``DL_POLY`` (or ``DLPOLY``\ ): this creates output in the ``TABLE`` format accepted by the ``DL_POLY`` simulation code.
* ``LAMMPS``: creates output suitable for use with the LAMMPS ``pair_style table`` (see :ref:`using-table-in-lammps`\ ).
* ``GULP`` produces a table for the ``GULP`` code by defining a set of separation, energy pairs using the ``GULP`` ``spline`` directive.

Defining table cut-off
++++++++++++++++++++++

The extent of the table is defined in the ``[Tabulation]`` section using the :ref:`nr <ref-potable-input-tabulation-nr>`\ , :ref:`dr <ref-potable-input-tabulation-dr>` and :ref:`cutoff <ref-potable-input-tabulation-cutoff>` options:

* ``dr`` defines the row increment (step-size) between table rows.
* ``cutoff`` gives the maximum separation to be tabulated.
* ``nr`` determines the number of rows in the tabulation.


Any two of ``nr``\ , ``dr``\ and ``cutoff`` can be used to define the extent and resolution of the tabulation. As an example all three of the following ``[Tabulation]`` sections would produce a table with 1000 rows, each separated by 0.01:

#. ``cutoff`` and ``nr``\ ::

	[Tabulation]
	target: LAMMPS
	cutoff : 9.99
	nr : 1000


#. ``cutoff`` and ``dr``\ ::

	[Tabulation]
	target: LAMMPS
	cutoff : 9.99
	dr : 0.01


#. ``nr`` and ``dr``\ ::

	[Tabulation]
	target: LAMMPS
	nr : 1000
	dr : 0.01
