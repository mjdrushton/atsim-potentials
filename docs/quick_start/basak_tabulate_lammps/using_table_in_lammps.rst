
.. _using-table-in-lammps:

********************************
Using ``Basak.lmptab`` in LAMMPS
********************************

`LAMMPS`_ input files are provided for use with the table file:

	* :download:`UO2.lmpstruct <UO2.lmpstruct>`: structure file for single UO\ :sub:`2` cell, that can be read with `read_data <http://lammps.sandia.gov/doc/read_data.html>`_ when `atom_style <http://lammps.sandia.gov/doc/atom_style.html>`_ ``full`` is used.
	* :download:`equilibrate.lmpin <equilibrate.lmpin>`: input file containing `LAMMPS`_ instructions. Performs 10ps of 300K NPT equilibration, creating  a 4×4×4 super-cell.

Copy these files into the same directory as ``Basak.lmptab``, the simulation can then be run using:

	.. code:: sh

		lammps -in equilibrate.lmpin -log equilibrate.lmpout -echo both 


The section of :download:`equilibrate.lmpin <equilibrate.lmpin>` which defines the potential model and makes use of the table file is as follows:

	.. literalinclude:: equilibrate.lmpin
		:lines: 14-26


	**Notes:**

		1. As `LAMMPS`_ uses ID numbers to define species the ``variable`` commands associate:
		   
		   * index 1 with variable ``$O`` 
		   * index 2 with ``$U`` to aid readability.

		3. The ``set type SPECIES_ID charge`` lines define the charges of oxygen and uranium.
		

		3. Uses the `hybrid/overlay`_ ``pair_style`` to combine the `coul/long`_ and `table`_ styles. 

			.. literalinclude:: equilibrate.lmpin
				:lines: 22

		   	* The `coul/long`_ style is used to calculate electrostatic interactions using the `pppm`_ ``kspace_style`` defined previously.
		   	* ``table linear 6500 pppm``: 

			   	* linear interpolation of table values should be used 
			   	* all 6500 rows of the table are employed
			   	* corrections appropriate to the `pppm`_ ``kspace_style`` will be applied.
		   

		4. 	Means that electrostatic interactions should be calculated between all pairs of ions.

		   	.. literalinclude:: equilibrate.lmpin
				:lines: 23
		  	
		
		5. Each ``pair_coeff`` reads an interaction from the ``Basak.lmptab`` file. 

		  	.. literalinclude:: equilibrate.lmpin
				:lines: 24-26
		   	
			* The general form is:

			   	* ``pair_coeff SPECIES_A SPECIES_A table TABLE_FILENAME TABLE_KEYWORD``
			   	* Here the  ``SPECIES_*`` use the ``$O`` and ``$U`` variables defined earlier.
			   	* ``TABLE_KEYWORD`` - the table file contains multiple blocks, each defining a single interaction. 
			   	* The ``TABLE_KEYWORD`` is the title of the block. The labels are of the form ``LABEL_A-LABEL_B`` albeit with the species sorted into alphabetical order.


.. _LAMMPS: http://lammps.sandia.gov
.. _hybrid/overlay: http://lammps.sandia.gov/doc/pair_hybrid.html
.. _coul/long: http://lammps.sandia.gov/doc/pair_coul.html
.. _table: http://lammps.sandia.gov/doc/pair_table.html
.. _pppm: http://lammps.sandia.gov/doc/kspace_style.html
