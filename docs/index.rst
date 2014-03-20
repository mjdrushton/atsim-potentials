.. atsim_potentials documentation master file, created by
   sphinx-quickstart on Thu Mar 13 22:54:46 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

****************
atsim.potentials
****************

Classical simulation codes typically contain a good selection of analytical forms for describing atomic interactions. Sometimes however, you may need to use a  potential that is not directly supported by the code. Luckily, most simulation codes allow you to provide tabulated potentials in which energies and forces, for a range of interatomic separations, are pre-calculated and specified as rows within a text file. The ``atsim.potentials`` package provides python modules to make the specification and tabulation of pair- and many-body potentials straightforward and consistent.

The following codes are supported for pair-potential tabulation:

* `LAMMPS`_ 
* `DL_POLY`_ 

    
Embedded Atom Model (EAM) potential tabulation is supported in the following formats:

* DYNAMO - as used by `LAMMPS`_ and several other codes.
* `DL_POLY`_


Contents
========

.. toctree::
    :maxdepth: 2

    installation
    quickstart
    potentials/pair_potential_tabulation
    potentials/eam_tabulation
    examples
    reference
    credits
    license

Contact
=======

``atsim.potentials`` was developed by Michael Rushton, if you have any problems, suggestions or queries please get in touch (contact details are available from my `homepage`_).

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. _LAMMPS: http://lammps.sandia.gov
.. _DL_POLY: http://www.stfc.ac.uk/cse/25526.aspx
.. _homepage: http://abulafia.mt.ic.ac.uk/groupmembers/michael