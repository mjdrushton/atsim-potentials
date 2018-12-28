***************************************************************************
`atsim.potentials` - Potential Model Tabulation for Atomic Scale Simulation
***************************************************************************

Classical simulation codes typically contain a good selection of analytical forms for describing atomic interactions. Sometimes however, you may need to use a  potential that is not directly supported by the code. Luckily, most simulation codes allow you to provide tabulated potentials in which energies and forces, for a range of interatomic separations, are pre-calculated and specified as rows within a text file. The ``atsim.potentials`` package provides python modules to make the specification and tabulation of pair- and many-body potentials straightforward and consistent.

Features
========

  * **Pair-Potential Tabulation:** Effective pair-potentials can be tabulated for multiple codes including:

    + `GULP`_
    + `LAMMPS`_ 
    + `DL_POLY`_ 

  * **Many-Bodied Potentials:** Embedded Atom Model (EAM) potential tabulation is supported in the following formats:
    + DYNAMO - as used by `LAMMPS`_ and several other codes:
   
      - Support for `LAMMPS`_ `eam`, `eam/fs`, `eam/alloy` pair-styles.

    + `DL_POLY`_\ : write `TABEAM` formatted files.

  * **No programming required:** ``atsim.potentials`` can be driven using its own potential definition format. Using simple configuration files complex models can be defined and tabulated without requiring any programming experience.
  * **Batteries included:** comes pre-loaded with a wide range of common-potential types.
  * **Potential Splining:** join different potentials together with splines.
  * **Flexible:** the ``atsim.potentials`` potential definition format allows the use of arbitrary mathematical formulae to define new potential functions. If this isn't sufficient it also provides a powerful Python API which should allow most tabulation tasks to be achieved.  


Contents
========

.. toctree::
    :maxdepth: 2

    installation
    quick_start/quickstart
    user_guide/user_guide
    reference/reference
    examples
    credits
    changes
    license

.. todolist::




Contact
=======

``atsim.potentials`` was developed by Michael Rushton, if you have any problems, suggestions or queries please get in touch at m.rushton@bangor.ac.uk

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. _LAMMPS: http://lammps.sandia.gov
.. _DL_POLY: https://www.scd.stfc.ac.uk/Pages/DL_POLY.aspx
.. _GULP: https://nanochemistry.curtin.edu.au/gulp/