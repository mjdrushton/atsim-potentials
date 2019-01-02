
.. _ref-potable-input-format:

************************
``potable`` input format
************************

.. contents::
	:local:

.. _ref-potable-input-tabulation:

[Tabulation]
============

.. _ref-potable-input-tabulation-cutoff:

cutoff
------

:Item: ``cutoff``
:Format: float
:Description: Defines upper bound of functions tabulated in terms of separation. This is used in a pair with ns tabulated in terms of separation. This directive is used together with :ref:`nr (number of rows) <ref-potable-input-tabulation-nr>` or :ref:`dr  (step size) <ref-potable-input-tabulation-dr>` to give the extent and resolution of a tabulated function.


.. _ref-potable-input-tabulation-dr:

dr
--

:Item: ``dr``
:Format: float
:Description: Defines the step size between rows of functions tabulated in terms of separation. This directive is used together with :ref:`nr <ref-potable-input-tabulation-nr>` or :ref:`cutoff <ref-potable-input-tabulation-cutoff>` to define resolution and extent of these functions.


.. _ref-potable-input-tabulation-nr:

nr
--

:Item: ``nr``
:Format: int
:Description: Defines the number of rows when functions are tabulated in terms of separation. This directive is used either with :ref:`dr <ref-potable-input-tabulation-dr>` or :ref:`cutoff <ref-potable-input-tabulation-cutoff>` to give the range and resolution of the tabulated function.

.. _ref-potable-input-tabulation-target:

target
------

:Item: ``target``
:Format: str
:Valid Options: ``DL_POLY|DLPOLY``, ``DLPOLY_EAM_fs`` ``DLPOLY_EAM``, ``GULP``, ``LAMMPS_eam_alloy|setfl``, ``LAMMPS``, ``setfl_fs``
:Description: Specifies the format that tabulation will be written in.