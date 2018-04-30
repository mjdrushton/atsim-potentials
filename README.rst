***********************************************************************************************
``atsim.potentials`` - python modules for working with potential models in atomistic simulation
***********************************************************************************************

Python modules and scripts to support atomistic simulation; ``atsim.potentials`` provides tools for working with pair and embedded atom method potential models. 

In particular includes functions for tabulation of pair and EAM models for `LAMMPS`_ and `DL_POLY`_.

Documentation
=============

Documentation, containing examples, API reference etc is hosted at http://atsimpotentials.readthedocs.org

Installation
============

Install Using Pip
-----------------

If you have `Pip <http://www.pip-installer.org/>`_ type the following to install ``atsim.potentials``:

.. code:: sh
	
	pip install atsim.potentials


Install from Source
-------------------

The source is hosted on `bitbucket`_ and can be cloned using `mercurial`_ as follows:

.. code:: sh

	hg clone https://bitbucket.org/mjdr/atsim_potentials  


alternatively a tarball of the source can be downloaded from the bitbucket downloads page `here <https://bitbucket.org/mjdr/atsim_potentials/downloads>`_ 

From the source directory install ``atsim.potentials`` using the following command:

.. code:: sh

	python setup.py install

Contact
=======

``atsim.potentials`` was developed by Michael Rushton, if you have any problems, suggestions or queries please get in touch at m.j.d.rushton@gmail.com .


License
=======

``atsim.potentials`` is licensed under the Apache 2.0 license. For more information,
please see LICENSE and NOTICE file.


.. _LAMMPS: http://lammps.sandia.gov
.. _DL_POLY: http://www.stfc.ac.uk/cse/25526.aspx
.. _bitbucket: http://https://bitbucket.org/mjdr/atsim_potentials/
.. _mercurial: http://mercurial.selenic.com
