
********************
Potential Tabulation
********************

Introduction
============

Classical simulation codes typically contain a good selection of analytical forms for describing atomic interactions. Sometimes however, you may need to use a potential form not directly supported by the code. Luckily, most simulation codes allow you to provide tabulated potentials in which energies and forces, for a range of interatomic separations, are pre-calculated and specified as rows within a text file. AtomsScripts provides python modules to make the specification and tabulation of pair- and many-body potentials straightforward.

The following codes are supported for pair-potential tabulation:

* `LAMMPS`_ 
* `DL_POLY`_ 
    

Embedded Atom Model (EAM) potential tabulation is supported in the following formats:

* DYNAMO - as used by `LAMMPS`_ and several other codes.
* `DL_POLY`_


Contents
========

.. toctree::
    :maxdepth: 3

    pair_potential_tabulation
    eam_tabulation



.. _LAMMPS: http://lammps.sandia.gov
.. _DL_POLY: http://www.stfc.ac.uk/cse/25526.aspx
