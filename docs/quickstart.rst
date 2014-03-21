.. _quick_start:

***********
Quick-Start
***********

The following example gives a complete python script showing how the potential API can be used to tabulate potentials for `DL_POLY`_ .

The following example (:download:`basak_tabulate.py`) shows how the UO\ :sub:`2` potential model of Basak [#basak03]_ can be tabulated:
    
    * the U + O interaction within this model combines Buckingham and Morse potential forms. Although `DL_POLY`_ natively supports both potential forms they cannot be combined with the code itself. By creating a ``TABLE`` file the Basak model can be described to `DL_POLY`_.
    * when executed from the command line this script will write tabulated potentials into a file named ``TABLE``.

.. literalinclude:: basak_tabulate.py


Tabulating the Potentials
-------------------------

Defining the Potentials
=======================

The first step to tabulating pair potentials is to define :class:`~atsim.potentials.Potential` objects (see :ref:`potential_objects`). Normally this involves creating a python function for the desired pair interaction before passing this to the :class:`atsim.potentials.Potential` together with labels name the species pair pertinent to the interaction. 

	* The functions ``f_OO`` and ``f_UU`` use the Buckingham form and are created using :func:`~atsim.potentials.potentialforms.buck` function factory (see :ref:`predefined_potential_forms` for more on the pre-defined forms provided):

		.. literalinclude:: basak_tabulate.py
			:lines: 7-15

	* The O-U interaction is a little more tricky to define as Buckingham and Morse potentials need to be combined. Pre-canned implementations of both of these are provided in :mod:`atsim.potentials.potentialforms` as :func:`~atsim.potentials.potentialforms.buck` and :func:`~atsim.potentials.potentialforms.morse`. Two functions are created, one for each component of the O-U interaction and stored in the ``buck_OU`` and ``morse_OU`` variables:

		.. literalinclude:: basak_tabulate.py
			:lines: 17-24

	* These are then composed into the desired function, ``f_OU``, using the :func:`~atsim.potentials.plus` function (see :ref:`combining_potential_forms`):

		.. literalinclude:: basak_tabulate.py
			:lines: 28


Make ``TABLE`` File
===================

The table file is written from the ``main()`` function of :download:`basak_tabulate.py`

	.. literalinclude:: basak_tabulate.py
		:pyobject: main


	* First the ``makePotentialObjects()`` function is called, returning a list of :class:`~atsim.potentials.Potential` that are stored in the ``potential_objects`` variable. 
	* The :func:`~atsim.potentials.writePotentials` function is then called with this list to create potentials with a maximum cut-off of 6.5Å and 6500 rows (i.e. a grid increment of 0.001 Å).

	* Now run the :download:`basak_tabulate.py` file (making sure you have :ref:`installed <installation>` ``atsim.potentials`` first):

		.. code::

			python basak_tabulate.py

	* This will create a DL_POLY ``TABLE`` file in the working directory. 


Using the ``TABLE`` File in DL_POLY
===================================

A set of DL_POLY files are provided allowing a simple NPT molecular dynamics equilibration simulation to be run against the ``TABLE`` file created in the previous step using :class:`~atsim.potentials.writePotentials`. Copy the files linked from the following table into the same directory as the ``TABLE`` file:

	======================================================		==============================================================
	File														Description
	------------------------------------------------------		--------------------------------------------------------------
	:download:`CONFIG <basak_tabulate_dlpoly/CONFIG>`			4×4×4 UO2:sub:`2` super-cell containing 768 atoms.
	:download:`CONTROL <basak_tabulate_dlpoly/CONTROL>`			Defines 300K equilibration run under NPT ensemble lasting 10ps.
	:download:`FIELD <basak_tabulate_dlpoly/FIELD>`				File defining potentials and charges.
	======================================================		==============================================================

* The ``FIELD`` file contains the directives relevant to the ``TABLE``  file:
  
	  .. literalinclude:: basak_tabulate_dlpoly/FIELD


* The following lines define the atom multiplicity and charges (O=-1.2\ *e* and U=2.4\ *e*):
	  
	  .. literalinclude:: basak_tabulate_dlpoly/FIELD
	  	 :lines: 5-9

* The ``vdw`` section states that the O-O, U-U and O-U interactions should be read from the ``TABLE`` file:

	  .. literalinclude:: basak_tabulate_dlpoly/FIELD
	  	 :lines: 10-14

* Once all the files are in the same directory, the simulation can be started by invoking ``DL_POLY``:
      	
    	.. code :: sh

    		DLPOLY.Z


.. _quick_start_lammps:

Quick-Start: LAMMPS
-------------------

Once the potential model has been defined as a series of :class:`~atsim.potentials.Potential` creating tabulations for different codes in different formats is fairly simple. The script described in this example is given in :download:`basak_tabulate_lammps.py`. This contains the same potential definition as the :ref:`previous example <quick_start>`, however the ``main()`` function has been modified to create a table suitable for `LAMMPS`_ :

.. literalinclude:: basak_tabulate_lammps.py
	:pyobject: main
	:emphasize-lines: 6,8

Only the two highlighted lines have been changed:

	1. the first changes the output filename to ``Basak.lmptab``
	2. the second has been changed from ``DL_POLY`` to ``LAMMPS`` in order to select the desired tabulation format.

Running the file creates the ``Basak.lmptab`` file:

	.. code:: sh

		python basak_tabulate_lammps.py


Using ``Basak.lmptab`` in LAMMPS
================================

`LAMMPS`_ input files are provided for use with the table file:

	* :download:`UO2.lmpstruct <basak_tabulate_lammps/UO2.lmpstruct>`: structure file for single UO:sub:`2` cell, that can be read with `read_data <http://lammps.sandia.gov/doc/read_data.html>`_ when `atom_style <http://lammps.sandia.gov/doc/atom_style.html>`_ ``full`` is used.
	* :download:`equilibrate.lmpin <basak_tabulate_lammps/equilibrate.lmpin>`: input file containing `LAMMPS`_ instructions. Performs 10ps of 300K NPT equilibration, creating  a 4×4×4 super-cell.

Copy these files into the same directory as ``Basak.lmptab``, the simulation can then be run using:

	.. code:: sh

		lammps -in equilibrate.lmpin -log equilibrate.lmpout -echo both 


The section of :download:`equilibrate.lmpin <basak_tabulate_lammps/equilibrate.lmpin>` which defines the potential model and makes use of the table file is as follows:

	.. literalinclude:: basak_tabulate_lammps/equilibrate.lmpin
		:lines: 14-26


	**Notes:**

		1. As `LAMMPS`_ uses ID numbers to define species the ``variable`` commands associate:
		   
		   * index 1 with variable ``$O`` 
		   * index 2 with ``$U`` to aid readability.

		3. The ``set type SPECIES_ID charge`` lines define the charges of oxygen and uranium.
		

		3. Uses the `hybrid/overlay`_ ``pair_style`` to combine the `coul/long`_ and `table`_ styles. 

			.. literalinclude:: basak_tabulate_lammps/equilibrate.lmpin
				:lines: 22

		   	* The `coul/long`_ style is used to calculate electrostatic interactions using the `pppm`_ ``kspace_style`` defined previously.
		   	* ``table linear 6500 pppm``: 

			   	* linear interpolation of table values should be used 
			   	* all 6500 rows of the table are employed
			   	* corrections appropriate to the `pppm`_ ``kspace_style`` will be applied.
		   

		4. 	Means that electrostatic interactions should be calculated between all pairs of ions.

		   	.. literalinclude:: basak_tabulate_lammps/equilibrate.lmpin
				:lines: 23
		  	
		
		5. Each ``pair_coeff`` reads an interaction from the ``Basak.lmptab`` file. 

		  	.. literalinclude:: basak_tabulate_lammps/equilibrate.lmpin
				:lines: 24-26
		   	
			* The general form is:

			   	* ``pair_coeff SPECIES_ID_1 SPECIES_ID_2 table TABLE_FILENAME TABLE_KEYWORD``
			   	* Here the  ``SPECIES_IDs`` use the ``$O`` and ``$U`` variables defines earlier.
			   	* ``TABLE_KEYWORD`` - the table file contains multiple blocks, each defining a single interaction. 
			   	* The ``TABLE_KEYWORD`` is the title of the block. The :func:`~atsim.potentials.writePotentials` function creates labels of the form ``LABEL_A-LABEL_B`` albeit with the species sorted into alphabetical order. This label format is described in greater detail :func:`here <atsim.potentials.writePotentials>`.



.. [#basak03] Basak, C. (2003). Classical molecular dynamics simulation of UO2 to predict thermophysical properties. *Journal of Alloys and Compounds*, **360** (1-2), 210–216. http://dx.doi.org/doi:10.1016/S0925-8388(03)00350-5
.. _DL_POLY: http://www.stfc.ac.uk/cse/25526.aspx
.. _LAMMPS: http://lammps.sandia.gov
.. _hybrid/overlay: http://lammps.sandia.gov/doc/pair_hybrid.html
.. _coul/long: http://lammps.sandia.gov/doc/pair_coul.html
.. _table: http://lammps.sandia.gov/doc/pair_table.html
.. _pppm: http://lammps.sandia.gov/doc/kspace_style.html
