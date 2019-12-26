.. _using-table-in-dlpoly:

***********************************
Using the ``TABLE`` File in DL_POLY
***********************************

A set of DL_POLY files are provided allowing a simple NPT molecular dynamics equilibration simulation to be run using a ``TABLE`` file created with ``atsim.potentials``. Copy the files linked from the following table into the same directory as the ``TABLE`` file:

	======================================================		==============================================================
	File														Description
	------------------------------------------------------		--------------------------------------------------------------
	:download:`CONFIG <../basak_tabulate_dlpoly/CONFIG>`			4×4×4 UO2:sub:`2` super-cell containing 768 atoms.
	:download:`CONTROL <../basak_tabulate_dlpoly/CONTROL>`			Defines 300K equilibration run under NPT ensemble lasting 10ps.
	:download:`FIELD <../basak_tabulate_dlpoly/FIELD>`				File defining potentials and charges.
	======================================================		==============================================================

* The ``FIELD`` file contains the directives relevant to the ``TABLE``  file:
  
	  .. literalinclude:: ../basak_tabulate_dlpoly/FIELD


* The following lines define the atom multiplicity and charges (O=-1.2\ *e* and U=2.4\ *e*):
	  
	  .. literalinclude:: ../basak_tabulate_dlpoly/FIELD
	  	 :lines: 5-9

* The ``vdw`` section states that the O-O, U-U and O-U interactions should be read from the ``TABLE`` file:

	  .. literalinclude:: ../basak_tabulate_dlpoly/FIELD
	  	 :lines: 10-14

* Once all the files are in the same directory, the simulation can be started by invoking ``DL_POLY``:
      	
    	.. code :: sh

    		DLPOLY.Z