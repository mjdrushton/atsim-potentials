.. _installation:

************
Installation
************

Install Using Pip
=================

If you have `Pip <http://www.pip-installer.org/>`_ type the following to install ``atsim.potentials``:

.. code:: sh
	
	pip install atsim.potentials


Install from Source
===================

The source is hosted on `bitbucket`_ and can be cloned using `mercurial`_ as follows:

.. code:: sh

	hg clone https://bitbucket.org/mjdr/atsim_potentials  


alternatively a tarball of the source can be downloaded from the bitbucket downloads page `here <https://bitbucket.org/mjdr/atsim_potentials/downloads>`_ 

From the source directory install ``atsim.potentials`` using the following command:

.. code:: sh

	python setup.py install


Build the Documentation
-----------------------

The documentation (which you are currently reading) can be built from source using (assuming sphinx is installed):

.. code:: sh

	python setup.py build_sphinx


This will place documents in ``build/sphinx/html`` within the source tree. 

Alternatively this documentation is hosted at http://atsimpotentials.readthedocs.org


.. _bitbucket: http://https://bitbucket.org/mjdr/atsim_potentials/
.. _mercurial: http://mercurial.selenic.com
