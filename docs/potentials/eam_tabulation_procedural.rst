
.. _python_api_procedural_eam_tabulation:

Procedural EAM Tabulation
=========================

As an alternative to the ``EAM_Tabulation`` objects described here :ref:`python_api_eam_tabulation` a procedural interface is also provided for EAM tabulation using the following functions:

===================================================  ===========  ===============  =================================================================
Function                                             File-Format  Simulation Code  Example
===================================================  ===========  ===============  =================================================================
:func:`~atsim.potentials.writeFuncFL`                ``funcfl``   `LAMMPS`_        :ref:`Example 1: Ag in LAMMPS <eam_example_1>` 
:func:`~atsim.potentials.writeSetFL`                 ``setfl``    `LAMMPS`_        :ref:`Example 2a: Al-Cu in LAMMPS <eam_example_2a>`
:func:`~atsim.potentials.writeTABEAM`                ``TABEAM``   `DL_POLY`_       :ref:`Example 2b: Al-Cu in LAMMPS <eam_example_2b>`
:func:`~atsim.potentials.writeSetFLFinnisSinclair`   ``setfl``    `LAMMPS`_        :ref:`Example 3a: Al-Fe Finnis-Sinclair in LAMMPS <eam_example_3a>`
:func:`~atsim.potentials.writeTABEAMFinnisSinclair`  ``TABEAM``   `DL_POLY`_       :ref:`Example 3b: Al-Fe Finnis-Sinclair in DL_POLY <eam_example_3b>`
===================================================  ===========  ===============  =================================================================


Examples
--------

.. _eam_example_1:

Example 1: Using :func:`~atsim.potentials.writeFuncFL` to Tabulate Ag Potential for LAMMPS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This example shows how to use :func:`~atsim.potentials.writeFuncFL` function to tabulate an EAM model for the simulation of Ag metal. How to use this tabulation within `LAMMPS`_ will then be demonstrated. The final tabulation script can be found in :download:`eam_tabulate_example1.py <eam_tabulate_example1.py>`.

The same model as used for the :class:`~atsim.potentials.eam_tabulation.SetFL_EAMTabulation` example (:ref:`eam_example_1_obj`) and is described the same way in python. In terms of the code, the only significant difference between the object based example and this one, is the use of :func:`~atsim.potentials.writeFuncFL` to tabulate the model into a file. The output format used in this example is also different, it uses the simpler ``funcfl`` format. Each ``funcfl`` file contains a single species, making alloy systems less convenient. Further more, alloy models are simulated by combining the ``funcfl`` files using pre-determined mixing rules meaning there is much less control over the specific interactions between the various elements in the alloy. To prodce the same ``setfl`` files as produced by the :class:`~atsim.potentials.eam_tabulation.SetFL_EAMTabulation` class, the :func:`~atsim.potentials.writeSetFL` function can be used (an example of which is given :ref:`below <eam_example_2a>`\ ).

The ``embed()`` and ``density()`` functions are defined for :math:`F_{\text{Ag}} (\rho)` and :math:`\rho_{\text{Ag}}` respectively:


.. literalinclude:: eam_tabulate_example1.py
  :lines: 2-13


The embedding and density functions should then be wrapped in an :class:`~atsim.potentials.EAMPotential` object to create a single item list:

.. literalinclude:: eam_tabulate_example1.py
  :lines: 23-24


Similarly the pair potential component, :math:`\phi_{\text{Ag}-\text{Ag}} (r_{ij}`, of the model can easily be defined as: 


.. literalinclude:: eam_tabulate_example1.py
  :pyobject: pair_AgAg
      

This can then be wrapped in a :class:`atsim.potentials.Potential` object to create a list of pair potentials. 

.. literalinclude:: eam_tabulate_example1.py
  :lines: 25


.. note:: 
  
  :func:`~atsim.potentials.writeFuncFL` only accepts a single :class:`~atsim.potentials.Potential` object and this should be the X-X interaction (where X is the species for which the ``funcfl`` tabulation is being created). Within the  'pair'-potential is tabulated as 

  :math:`\sqrt{\frac{\phi(r_{ij}) r_{ij}}{27.2 \times 0.529}}`

  The numerical constants (27.2 and 0.529) convert from eV and Å into the units of Hartree and Bohr radius used by the ``funcfl`` format. The square rooting of the potential function is important: the simulation code effectively reconstitutes a pair potential by multiplying two of these tabulated square-rooted functions (one for each species in each interacting pair) together. If atoms :math:`i` and :math:`j` in an interacting pair, have the same species then effectively the original pair-potential is obtained (albeit multiplied by :math:`r_{ij}`).

  By comparison, if multiple ``funcfl`` files are used to define multiple species within a simulation (e.g. for alloy systems), then the pair potential functions of each species are effectively 'mixed' when they are multiplied together for heterogeneous atom pairs. If more control is required, with pair-potential functions specific to distinct pairs of species being necessary, then the ``setfl`` format produced by the :func:`~atsim.potentials.writeSetFL` and :func:`~atsim.potentials.writeSetFLFinnisSinclair` functions should be used instead.


Now all the components of the model have been defined a table file can be created in the ``funcfl`` format. Before doing this, it is necessary to choose  appropriate density and separation cut-offs together with :math:`d r_{ij}` and :math:`d \rho` increments for the density/pair functions and embedding function respectively:

  * Here a :math:`d \rho` value of 0.001 will be used and 50000 density values tabulated.
  * This means the maximum density that can be accepted by the embedding function is :math:`49999 \times 0.001 = 49.999`
  * :math:`dr = 0.001` Å using 12000 rows.
  * The pair-potential cut-off and the maximum :math:`r_{ij}` value for the density function is therefore 11.999 Å.

Invoking the :func:`~atsim.potentials.writeFuncFL` function with these values and the :class:`~atsim.potentials.EAMPotential` and :class:`.~atsim.potentialsPotential` objects, can be used to tabulate the Ag potential into the ``Ag.eam`` file:
    
.. literalinclude:: eam_tabulate_example1.py
  :language: python
  :lines: 27-42


Putting this together the following script is obtained (this script can also be downloaded :download:`eam_tabulate_example1.py <eam_tabulate_example1.py>`:

.. literalinclude:: eam_tabulate_example1.py
  :language: python

  
Running this script will produce a table file named ``Ag.eam`` in the same directory as the script:

.. code:: sh

  python eam_tabulate_example1.py


Using the ``Ag.eam`` file within LAMMPS
"""""""""""""""""""""""""""""""""""""""

This section of the example will now demonstrate how the table file can be used used to perform a static energy minimisation of an FCC Ag structure in LAMMPS.

Place the following in a file called :download:`fcc.lmpstruct <fcc.lmpstruct>` in the same directory as the ``Ag.eam`` file you created previously. This describes a single FCC cell with a wildly inaccurate lattice parameter:

.. literalinclude:: fcc.lmpstruct


The following LAMMPS input file describes a minimisation run. The lines describing potentials are highlighted. Put its contents in a file called :download:`example1_minimize.lmpin <example1_minimize.lmpin>`:

.. literalinclude:: example1_minimize.lmpin
  :emphasize-lines: 7,8


The ``pair_style eam`` command tells LAMMPS to use the EAM and expect ``pair_coeff`` commands mapping atom types to particular table files:

.. literalinclude:: example1_minimize.lmpin
  :lines: 7
  
The following ``pair_coeff`` directive indicates that the interaction between atom-type 1 (Ag) with itself should use the ``funcfl`` formatted file contained within ``Ag.eam``:

.. literalinclude:: example1_minimize.lmpin
  :lines: 8

The example can then be run by invoking LAMMPS:

.. code:: sh

  lammps -in example1_minimize.lmpin



.. _eam_example_2a:

Example 2a: Tabulate Al-Cu Alloy Potentials Using :func:`~atsim.potentials.writeSetFL` for LAMMPS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Within the following example the process required to generate and use a ``setfl`` file that tabulates the Al-Cu alloy model of Zhou et al [#Zhou2004]_.  By comparison to the ``funcfl`` format, ``setfl`` allows multiple elements to be given in the same file and additionally pair-potentials for particular pairs of interacting species can be specified (``funcfl`` relies on the simulation code to 'mix' pair-potentials within alloy systems). The :download:`eam_tabulate_example2a.py` gives a complete example of how the Zhou model can be tabulated.

This example is almost entirely the same as that given for the object based interface (:ref:`eam_example_2a_obj`) with the only difference being the use of the :func:`~atsim.potentials.writeSetFL` function for the final tabulation. For a description of the Zhou model and how it is coded in python please :ref:`see here <eam_example_2a_obj_model_description>`\ .

Putting everything together gives the following script (which can also be downloaded using the following link :download:`eam_tabulate_example2a.py <eam_tabulate_example2a.py>`:). Running this (``python eam_tabulate_example2a.py``) produces the ``Zhou_AlCu.eam.alloy`` file in current working directory.

.. literalinclude:: eam_tabulate_example2a.py
  :language: python


.. seealso::

  * See :ref:`eam_example_2a_obj_using_zhou` for details of how to use the tabulation file with `LAMMPS`_\ .


.. _eam_example_2b:

Example 2b: Tabulate Al-Cu Alloy Potentials Using :func:`~atsim.potentials.writeTABEAM` for DL_POLY
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The tabulation script used with :ref:`Example 2a <eam_example_2a>` can be easily modified to produce the ``TABEAM`` format expected by the `DL_POLY`_ simulation code by using the :func:`~atsim.potentials.writeTABEAM`\ . See the tabulation script for this example: :download:`eam_tabulate_example2b.py`.

.. literalinclude:: eam_tabulate_example2b.py
  :pyobject: main


.. seealso::

  * See the object oriented version of this example :ref:`eam_example_2b_obj`\ .




.. _eam_example_3a:

Example 3a: Tabulate Al-Fe Finnis-Sinclair Potentials Using :func:`~atsim.potentials.writeSetFLFinnisSinclair` for LAMMPS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This example will show how to reproduce the EAM model described by Mendelev et al. for Fe segregation at grain boundaries within Al [#Mendelev2011]_. As a result this example effectively shows how to reproduce the ``AlFe_mm.eam.fs`` file provided with the LAMMPS source distribution using the :func:`~atsim.potentials.writeSetFLFinnisSinclair` function.

The example uses the :func:`~atsim.potentials.writeSetFLFinnisSinclair` function to produce files supported by the LAMMPS ``pair_style eam/fs`` command. 

The potential model and definition of potential objects is detailed in :ref:`eam_example_3b_obj` which uses a tabulation class but is otherwise very similar to this example. Having defined the list of :class:`~atsim.potentials.EAMPotential` instances the :func:`~atsim.potentials.writeSetFLFinnisSinclair` function is called, in this case writing the data to ``Mendelev_Al_Fe.eam.fs`` in the current directory:

.. literalinclude:: eam_tabulate_example3a.py
  :language: python
  :pyobject: main


The full tabulation script can be downloaded as :download:`eam_tabulate_example3a.py`\ .

.. _eam_example_3b:

Example 3b: Tabulate Al-Fe Finnis-Sinclair Potentials Using :func:`~atsim.potentials.writeTABEAMFinnisSinclair` for DL_POLY
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Using exactly the same model definition as for :ref:`Example 3a <eam_example_3a>`, the Al-Fe model can be re-tabulated for DL_POLY with minimal modification to the ``main()`` function. The modified version of the tabulation script can be found in :download:`eam_tabulate_example3b.py`.  

The ``main()`` function is given below:

.. literalinclude:: eam_tabulate_example3b.py
  :pyobject: main

Excluding the import statement at the top of the file, only two lines have been changed (highlighted). The first changes the filename to ``TABEAM`` whilst the second tells python to call :func:`~atsim.potentials.writeTABEAMFinnisSinclair` instead of :func:`~atsim.potentials.writeSetFLFinnisSinclair`\ :

.. literalinclude:: eam_tabulate_example3b.py
  :lines: 129-136
  :emphasize-lines: 2



.. rubric:: Footnotes

.. [#sutton1990] A.P. Sutton,  and  J. Chen,  "Long-range Finnis-Sinclair potentials", *Philos. Mag. Lett.* **61** (1990) 139 `doi:10.1080/09500839008206493 <https://dx.doi.org/10.1080/09500839008206493>`_\ .
.. [#Zhou2004] X. Zhou, R. Johnson and H. Wadley, "Misfit-energy-increasing dislocations in vapor-deposited CoFe/NiFe multilayers", *Phys. Rev. B.* **69** (2004) 144113.  
.. [#Mendelev2011] M.I. Mendelev, D.J. Srolovitz, G.J. Ackland, and S. Han, "Effect of Fe Segregation on the Migration of a Non-Symmetric Σ5 Tilt Grain Boundary in Al", *J. Mater. Res.* **20** (2011) 208.


.. _DL_POLY: http://www.stfc.ac.uk/cse/25526.aspx
.. _LAMMPS: http://lammps.sandia.gov