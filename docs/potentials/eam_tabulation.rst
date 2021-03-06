
.. toctree::
  :hidden:

  eam_tabulation_procedural.rst



.. _python_api_eam_tabulation:

Embedded Atom Method (EAM) Tabulation
=====================================

An EAM model is defined by constructing instances of  :class:`atsim.potentials.EAMPotential` describing each species within the model. :class:`~atsim.potenials.EAMPotential` encapsulates the density and embedding functions specific to each species' many bodied interactions. In addition the purely pairwise interactions within the EAM are defined using a list of :class:`atsim.potentials.Potential` objects.


Once the EAM model has been described in terms of  :class:`EAMPotential <atsim.potentials.EAMPotential>` and :class:`Potential <atsim.potentials.Potential>` objects it can be tabulated for specific simulation codes. This is done by using the ``EAM_Tabulation`` objects from the :mod:`atsim.potentials.eam_tabulation` module:

+-------------------------------------------------------------------------------+------------------------------------------------------------------------+-----------------+---------------------------------------+
| Class                                                                         | Format                                                                 | Simulation Code | Example                               |
+===============================================================================+========================================================================+=================+=======================================+
| :class:`~atsim.potentials.eam_tabulation.SetFL_EAMTabulation`                 | | ``setfl``\ ,                                                         | LAMMPS          | | :ref:`eam_example_1_obj`            |
|                                                                               | | `pair_style eam/alloy <http://lammps.sandia.gov/doc/pair_eam.html>`_ |                 | |                                     |
|                                                                               |                                                                        |                 | | :ref:`eam_example_2a_obj`           |
+-------------------------------------------------------------------------------+------------------------------------------------------------------------+-----------------+---------------------------------------+
| :class:`~atsim.potentials.eam_tabulation.SetFL_FS_EAMTabulation`              | `pair_style eam/fs <http://lammps.sandia.gov/doc/pair_eam.html>`_      | LAMMPS          | :ref:`eam_example_3a_obj`             |
+-------------------------------------------------------------------------------+------------------------------------------------------------------------+-----------------+---------------------------------------+
| :class:`~atsim.potentials.eam_tabulation.TABEAM_EAMTabulation`                | ``TABEAM``                                                             | DL_POLY         | :ref:`eam_example_2b_obj`             |
+-------------------------------------------------------------------------------+------------------------------------------------------------------------+-----------------+---------------------------------------+
| :class:`~atsim.potentials.eam_tabulation.TABEAM_FinnisSinclair_EAMTabulation` | ``EEAM TABEAM``                                                        | DL_POLY         | :ref:`eam_example_3b_obj`             |
+-------------------------------------------------------------------------------+------------------------------------------------------------------------+-----------------+---------------------------------------+
| :class:`~atsim.potentials.eam_tabulation.Excel_EAMTabulation`                 | ``.xlsx``                                                              |                 |                                       |
+-------------------------------------------------------------------------------+------------------------------------------------------------------------+-----------------+---------------------------------------+
| :class:`~atsim.potentials.eam_tabulation.Excel_FinnisSinclair_EAMTabulation`  | ``.xlsx``                                                              |                 |                                       |  
+-------------------------------------------------------------------------------+------------------------------------------------------------------------+-----------------+---------------------------------------+


Even though the use of ``EAM_Tabulation`` objects is preferred a legacy procedural interface is also provided. This is described here: :ref:`python_api_procedural_eam_tabulation`\ .


Examples
--------

.. _eam_example_1_obj:

Example 1: Using :class:`~atsim.potentials.eam_tabulation.SetFL_EAMTabulation` to Tabulate Ag Potential for LAMMPS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This example shows how to use the :class:`~atsim.potentials.eam_tabulation.SetFL_EAMTabulation` class to tabulate an EAM model for the simulation of Ag metal. How to use this tabulation within `LAMMPS`_ will then be demonstrated. The final tabulation script can be found in :download:`eam_example1.py <eam_example1.py>`.

.. seealso::

  * A ``potable`` version of this example is given here: :ref:`potable_sutton_ag_example`\ .

.. _eam_example_1_obj_model_description:

Model Description
"""""""""""""""""

Within this example the Ag potential of Sutton will be tabulated [#sutton1990]_. Within the EAM the energy (:math:`E_i`) of an atom :math:`i` whose species is :math:`\alpha` is given by:

.. math::
  E_i = F_\alpha \left( \sum_{j \neq i} \rho_\beta(r_{ij}) \right) + \frac{1}{2} \sum_{j \neq i} \phi_{\alpha \beta} (r_{ij})


.. note::

  * :math:`\rho_\beta(r_{ij})` is the density function which gives the electron density for atom :math:`j` with species :math:`\beta` as a function of its separation from atom :math:`i`, :math:`r_{ij}`.
  * The electron density for atom :math:`i` is obtained by summing over the density (:math:`\rho_\beta (r_{ij}`) contributions due to its neighbours.
  * The embedding function :math:`F_\alpha(\rho)` is used to calculate the many-bodied energy contribution from this summed electron density.
  * The sum :math:`\frac{1}{2} \sum_{j \neq i} \phi_{\alpha \beta} (r_{ij})` gives the pair-potential contribution to atom :math:`i`'s energy. 
  * :math:`\phi_{\alpha \beta} (r_{ij})` are simply pair potentials that describe the energy between two atoms as a function of their separation.


The embedding function used by Sutton is:

.. math::
  F_\alpha (\rho)  = - \sqrt{\rho}

and the density function is:

.. math::
  \rho_\beta (r_{ij}) = \left( \frac{a}{r_{ij}} \right)^m

whilst pair interactions are given by:

.. math::
  \phi_{\alpha \beta} (r_{ij}) = \left( \frac{b}{r_{ij}} \right)^n
  

The model parameters are given as:

============ =====================================================
Parameter    Value
============ =====================================================
:math:`m`    6
:math:`n`    12
:math:`a`    :math:`2.928323832 \text{Å} \text{eV}^{\frac{1}{3}}`
:math:`b`    :math:`2.485883762 \text{eV}^\frac{1}{12} \text{Å}`
============ =====================================================

Define the Model
""""""""""""""""

It is now necessary to describe the model in python code. Hard-coding the model parameters from the previous table, ``embed()`` and ``density()`` functions can be defined for :math:`F_{\text{Ag}} (\rho)` and :math:`\rho_{\text{Ag}}` respectively:


.. literalinclude:: eam_example1.py
  :lines: 2-16


The embedding and density functions should then be wrapped in an :class:`~atsim.potentials.EAMPotential` object to create a single item list:

.. literalinclude:: eam_example1.py
  :lines: 26


Similarly the pair potential component, :math:`\phi_{\text{Ag}-\text{Ag}} (r_{ij})`, of the model can easily be defined as: 


.. literalinclude:: eam_example1.py
  :pyobject: pair_AgAg
      

This can then be wrapped in a :class:`atsim.potentials.Potential` object to create a list of pair potentials. 

.. literalinclude:: eam_example1.py
  :lines: 27



Now all the components of the model have been defined a table file can be created in the ``setfl`` format. Before doing this, it is necessary to choose  appropriate density and separation cut-offs together with the number of rows in the density/pair functions (``nr``) and embedding function (``nrho``) respectively:

  * Here 50000 density values will be tabulated to a cutoff of 50.0.
  * The pair-potential cut-off and the maximum :math:`r_{ij}` value for the density function is 12 Å and both will have 12000 rows.

An instance of :class:`atsim.potentials.eam_tabulation.SetFL_EAMTabulation` is created with the :class:`~atsim.potentials.EAMPotential` and :class:`~atsim.potentials.Potential` objects. This object is then used to tabulate the Ag potential by calling the :meth:`~atsim.potentials.eam_tabulation.SetFL_EAMTabulation.write` method with the ``Ag.eam.alloy`` file object:
    
.. literalinclude:: eam_example1.py
  :language: python
  :lines: 27-42


Putting this together the following script is obtained (this script can also be downloaded :download:`eam_example1.py <eam_example1.py>`:

.. literalinclude:: eam_example1.py
  :language: python

  
Running this script will produce a table file named ``Ag.eam.alloy`` in the same directory as the script:

.. code:: sh

  python eam_example1.py


Using the ``Ag.eam.alloy`` file within LAMMPS
"""""""""""""""""""""""""""""""""""""""""""""

This section of the example will now demonstrate how the table file can be used used to perform a static energy minimisation of an FCC Ag structure in LAMMPS.

Place the following in a file called :download:`fcc.lmpstruct <fcc.lmpstruct>` in the same directory as the ``Ag.eam.alloy`` file you created previously. This describes a single FCC cell with a wildly inaccurate lattice parameter:

.. literalinclude:: fcc.lmpstruct


The following LAMMPS input file describes a minimisation run. The lines describing potentials are highlighted. Put its contents in a file called :download:`example_eam_alloy_minimize.lmpin <example_eam_alloy_minimize.lmpin>`:

.. literalinclude:: example_eam_alloy_minimize.lmpin
  :emphasize-lines: 7,8


The ``pair_style eam/alloy`` command tells LAMMPS to use the EAM and expect ``pair_coeff`` commands mapping atom types to particular table files:

.. literalinclude:: example_eam_alloy_minimize.lmpin
  :lines: 7
  
The following ``pair_coeff`` directive indicates that the interaction between atom-type 1 (Ag) with itself should use the ``setfl`` formatted file contained within ``Ag.eam.alloy``\ . The ``Ag`` label at the end of the line indicates that atom-type 1 should be associated with this label in the table file:

.. literalinclude:: example_eam_alloy_minimize.lmpin
  :lines: 8

The example can then be run by invoking LAMMPS:

.. code:: sh

  lammps -in example_eam_alloy_minimize.lmpin



.. _eam_example_2a_obj:

Example 2a: Tabulate Al-Cu Alloy Potentials Using :class:`~atsim.potentials.eam_tabulation.SetFL_EAMTabulation` for LAMMPS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Within the following example the process required to generate and use a ``setfl`` file that tabulates the Al-Cu alloy model of Zhou et al [#Zhou2004]_.  In comparison to the previous example this example contains density and embedding functions for multiple elements and includes pair-potentials specific to pairs of interacting species. The :download:`eam_example2a.py` file gives a complete example of how the Zhou model can be tabulated.

.. _eam_example_2a_obj_model_description:

Model Description
"""""""""""""""""

The model makes use of the EAM as described above (see Example 1 :ref:`eam_example_1_obj_model_description`) . The density function, :math:`\rho_\beta (r_{ij})` for an atom :math:`j` of species :math:`\beta` separated from atom :math:`i` by :math:`r_{ij}` is:

.. math::
  
  \rho_\beta (r_{ij}) = \frac{f_e \exp \left[ -\omega (r_{ij} / r_e - 1) \right] }{1+(r_{ij}/r_e - \lambda)^{20}}

where :math:`f_e`, :math:`r_e`, :math:`\omega` and :math:`\lambda` are parameters specific to species :math:`\beta`. The pair-potential function acting between species :math:`\alpha`--:math:`\beta` is obtained by combining the density functions of the interacting species:

.. math::

  \phi_{\alpha\beta} (r_{ij}) = \frac{1}{2} \left[ \frac{\rho_\beta(r_{ij})}{\rho_\alpha(r_{ij})} \phi_{\alpha\alpha}(r_{ij}) + \frac{\rho_\alpha(r_{ij})}{\rho_\beta(r_{ij})} \phi_{\beta\beta}(r_{ij}) \right]

The homogeneous pair-interactions, :math:`\phi_{\alpha\alpha}(r_{ij})` and :math:`\phi_{\beta\beta}(r_{ij})` have the form:

.. math::

  \phi_{\alpha \alpha}(r_{ij}) = \frac{A \exp \left[ -\gamma (r_{ij} / r_e - 1) \right] }{1+(r_{ij}/r_e - \kappa)^{20}} - \frac{B \exp \left[ -\omega (r_{ij} / r_e - 1) \right] }{1+(r_{ij}/r_e - \lambda)^{20}}

again, :math:`A`, :math:`B`, :math:`\gamma`, :math:`\omega`, :math:`\kappa` and :math:`\omega` are parameters specific to the species :math:`\alpha`.

The embedding function for each species, :math:`F_\alpha (\rho)`, is defined over three density ranges using the following:

.. math::

  F_\alpha (\rho) = \left\{ 
    \begin{array}{ll}
      \sum_{i=0}^3 F_{ni}\left( \frac{\rho}{\rho_n} - 1\right)^i & \rho < \rho_n, & \rho_n = 0.85 \rho_e \\
      \sum_{i=0}^3 F_{i}\left( \frac{\rho}{\rho_e} - 1\right)^i & \rho_n \leq \rho < \rho_0, & \rho_0 = 1.15 \rho_e \\
      F_e \left[1-\eta\ln\left(\frac{\rho}{\rho_s}\right)\right]\left(\frac{\rho}{\rho_s}\right)^\eta &\rho_0 \leq \rho & \\
    \end{array}
  \right.

The model parameters for Cu and Al are given in the following table:

    ================ ========== =========
    Parameter        Cu         Al
    ================ ========== =========
    :math:`r_e`      2.556162   2.863924
    :math:`f_e`      1.554485   1.403115
    :math:`\rho_e`   21.175871  20.418205
    :math:`\rho_s`   21.175395  23.195740
    :math:`\gamma`   8.127620   6.613165
    :math:`\omega`   4.334731   3.527021
    :math:`A`        0.396620   0.314873
    :math:`B`        0.548085   0.365551
    :math:`\kappa`   0.308782   0.379846
    :math:`\lambda`  0.756515   0.759692
    :math:`F_{n0}`   -2.170269  -2.807602
    :math:`F_{n1}`   -0.263788  -0.301435
    :math:`F_{n2}`   1.088878   1.258562
    :math:`F_{n3}`   -0.817603  -1.247604
    :math:`F_0`      -2.19      -2.83
    :math:`F_1`      0          0
    :math:`F_2`      0.561830   0.622245
    :math:`F_3`      -2.100595  -2.488244
    :math:`\eta`     0.310490   0.785902
    :math:`F_e`      -2.186568  -2.824528
    ================ ========== =========

.. note:: 
  The Al :math:`A` value is given as 0.134873 in Zhou's original *Phys. Rev. B* paper. However parameter file provided by Zhou for this model, at http://www.ctcms.nist.gov/potentials/Zhou04.html gives the parameter as 0.314873. It is this latter value that is used here.

  In addition the final term of the embedding function has been modified to match that used in fortran tabulation code also provided at http://www.ctcms.nist.gov/potentials/Zhou04.html


Define the Model
""""""""""""""""

A series of python functions are defined to describe the embedding, density and pair interaction functions. To encourage code re-use a number of function factories are defined. Using the parameters passed to them they return specialised functions appropriate for the parameters. The given factory functions make use of python's support for `closures <http://www.shutupandship.com/2012/01/python-closures-explained.html>`_  in their implementation.

The ``makeFunc()`` factory function is used to define density functions. As this functional form is also used as a component of the pair-potentials ``makeFunc()`` is re-used within the ``makePairPotAA()`` factory function.

.. literalinclude:: eam_example2a.py
  :pyobject: makeFunc

The following factory returns the functions used to describe the homogeneous Al-Al and Cu-Cu pair-potential interactions:

.. literalinclude:: eam_example2a.py
  :pyobject: makePairPotAA


Whilst ``makePairPotAB()`` describes the Al-Cu pair-potential:

.. literalinclude:: eam_example2a.py
  :pyobject: makePairPotAB


The ``makeEmbed()`` function describes the embedding function:

.. literalinclude:: eam_example2a.py
  :pyobject: makeEmbed


Lists of :class:`~atsim.potentials.EAMPotential` and :class:`~atsim.potentials.Potential` objects are created and returned as a tuple by the ``makePotentialObjects()`` function within :download:`eam_example2a.py`. Before invoking the factory functions we just defined, the model parameters are assigned to easily identifiable variables within this function:

.. literalinclude:: eam_example2a.py
  :lines: 58-121


Now the functions required by the :class:`~atsim.potentials.EAMPotential` instances for Al and Cu can be created:

.. literalinclude:: eam_example2a.py
  :lines: 93-99


Now these are wrapped up in :class:`~atsim.potentials.EAMPotential` objects to give the ``eamPotentials`` list:

.. literalinclude:: eam_example2a.py
  :lines: 102-104


Similarly, using the ``makePairPotAA()`` and ``makePairPotAB()`` function factories the :class:`~atsim.potentials.Potential` objects required for the tabulation are defined:

.. literalinclude:: eam_example2a.py
  :lines: 106-119


Now we have all the objects required for  :class:`~atsim.potentials.eam_tabulation.SetFL_EAMTabulation`\ . The next excerpt calls ``makeObjects()`` to get the EAM and pair-potential objects before creating the tabulation object, and invoking its :meth:`~atsim.potentials.eam_tabulation.SetFL_EAMTabulation.write` method to write the data into a file called ``Zhou_AlCu.eam.alloy``:

.. literalinclude:: eam_example2a.py
  :pyobject: main


Putting this all together gives the following script (which can also be downloaded using the following link :download:`eam_example2a.py <eam_example2a.py>`:). Running this (``python eam_example2a.py``) produces the ``Zhou_AlCu.eam.alloy`` file in current working directory.

.. literalinclude:: eam_example2a.py
  :language: python


.. _eam_example_2a_obj_using_zhou:

Using the ``Zhou_AlCu.eam.alloy`` file within LAMMPS
""""""""""""""""""""""""""""""""""""""""""""""""""""

Within LAMMPS the ``setfl`` files generated by :class:`~atsim.potentials.eam_tabulation.SetFL_EAMTabulation` are used with the `eam/alloy <http://lammps.sandia.gov/doc/pair_eam.html>`_ pair_style. The ``pair_coeff`` directive used with this ``pair_style`` effectively maps LAMMPS species numbers to the element names within the table file.

Single Element Systems 
:::::::::::::::::::::::

Assuming a LAMMPS system containing only Al (i.e. Al is species 1) then the ``pair_style`` and ``pair_coeff`` directives would be given as:

::

  pair_style eam/alloy
  pair_coeff * * Zhou_AlCu.eam.alloy Al

Likewise if a copper system was being simulated:

::

  pair_style eam/alloy
  pair_coeff * * Zhou_AlCu.eam.alloy Cu


Mixed Al-Cu System
::::::::::::::::::

For an Al-Cu system where Al is species 1 and Cu species 2 then the directives would be:

::

  pair_style eam/alloy
  pair_coeff * * Zhou_AlCu.eam.alloy Al Cu


Or if Cu was 1 and Al 2:

::

  pair_style eam/alloy
  pair_coeff * * Zhou_AlCu.eam.alloy Cu Al


.. _eam_example_2b_obj:

Example 2b: Tabulate Al-Cu Alloy Potentials Using :class:`~atsim.potentials.eam_tabulation.TABEAM_EAMTabulation` for DL_POLY
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The tabulation script used with :ref:`Example 2a <eam_example_2a_obj>` can be easily modified to produce the ``TABEAM`` format expected by the `DL_POLY`_ simulation code. See the tabulation script for this example: :download:`eam_example2b.py`.

The :class:`~atsim.potentials.EAMPotential` and :class:`~atsim.potentials.Potential` lists are created in exactly the same way as :ref:`Example 2a <eam_example_2a>`, however rather than creating an instance of :class:`~atsim.potentials.eam_tabulation.SetFL_EAMTabulation` in the ``main()`` function it is modified to use the DL_POLY specific :class:`~atsim.potentials.eam_tabulation.TABEAM_EAMTabulation` class instead and to write into a file named ``TABEAM``. The ``main()`` function of :download:`eam_example2b.py` is now given:

.. literalinclude:: eam_example2b.py
  :pyobject: main


Using the ``TABEAM`` file with `DL_POLY`_
"""""""""""""""""""""""""""""""""""""""""

Running :download:`eam_example2b.py` will create a file named ``TABEAM`` in the working directory. This should be copied into the simulation directory containing the `DL_POLY`_ input files (``CONTROL``, ``CONFIG`` and ``FIELD``). 

The following should be added at the bottom of the ``FIELD`` file:

::

  metal 3
  Al Al eam
  Cu Cu eam
  Al Cu eam



.. _eam_example_3a_obj:

Example 3a: Tabulate Al-Fe Finnis-Sinclair Potentials Using :class:`~atsim.potentials.eam_tabulation.SetFL_FS_EAMTabulation` for LAMMPS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This example will show how to reproduce the EAM model described by Mendelev et al. for Fe segregation at grain boundaries within Al [#Mendelev2011]_. As a result this example effectively shows how to reproduce the ``AlFe_mm.eam.fs`` file provided with the LAMMPS source distribution using the :class:`~atsim.potentials.eam_tabulation.SetFL_FS_EAMTabulation` class.

.. seealso::

  * :ref:`Finnis-Sinclair example using potable <many_body_models_potable_finnis_sinclair>`
  * :ref:`eam_example_3a`


The file format created by :class:`atsim.potentials.eam_tabulation.SetFL_FS_EAMTabulation` is supported by the LAMMPS ``pair_style eam/fs`` command. This adds an additional level of flexibility in comparison to the ``eam/alloy`` style; when calculating the density surrounding an atom with species :math:`\alpha`, each neighbouring atom's contribution to the density is calculated as a function of its separation from the central atom using :math:`\rho_{\alpha\beta}(r_{ij})`. This means that the density function is now specific to both the central atom species, :math:`\alpha` **and** that of the surrounding atom, :math:`\beta`. By comparison when using ``eam/alloy`` tabulations the same :math:`\rho_\beta(r_{ij})` function is used, no matter the type of the central atom. This means that the equation describing ``eam/fs`` style models becomes:


.. math::
  E_i = F_\alpha \left( \sum_{j \neq i} \rho_{\alpha\beta}(r_{ij}) \right) + \frac{1}{2} \sum_{j \neq i} \phi_{\alpha \beta} (r_{ij})


Here a binary Al, Fe,  model is being described and the resultant ``eam/fs`` file should contain definitions for the following:

* **Pair-Potentials**: :math:`\phi_{\text{Al}\text{Al}}(r_{ij})`, :math:`\phi_{\text{Fe}\text{Fe}}(r_{ij})` and :math:`\phi_{\text{Al}\text{Fe}}(r_{ij})`.
* **Embedding-Functions**: :math:`F_{\text{Al}}(\rho)` and :math:`F_{\text{Fe}}(\rho)`. 
* **Density-Functions**: :math:`\rho_{\text{Al}\text{Al}}(r_{ij})`, :math:`\rho_{\text{Fe}\text{Fe}}(r_{ij})`, :math:`\rho_{\text{Al}\text{Fe}}(r_{ij})` and :math:`\rho_{\text{Fe}\text{Al}}(r_{ij})`.
  
From this it can be seen that, when using ``eam/fs`` style potentials, the density functions must have both the :math:`\alpha\beta` **and** :math:`\beta\alpha` interactions specified.

Although both the :math:`\alpha\beta` and :math:`\beta\alpha` can be described using ``eam/fs`` files, the Mendelev model used in this example uses the same density function for both Al-Fe and Fe-Al cross density functions [#Mendelev2011]_. 

Using :class:`~atsim.potentials.eam_tabulation.SetFL_FS_EAMTabulation` to Tabulate the Model
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

As in previous examples it is necessary to define pair, density and embedding functions in python code that are then wrapped in :class:`~atsim.potentials.EAMPotential` and :class:`~atsim.potentials.Potential` objects to be passed to the tabulation function. For brevity only the names of the functions, as defined in the attached example file (:download:`eam_example3a.py`) are now  given:

* **Pair-Potentials:**

  * ``def ppfuncAlAl(r):`` - Al-Al pair-potential :math:`\phi_{\text{Al}\text{Al}}(r_{ij})`.
  * ``def ppfuncAlFe(r):`` - Al-Fe pair-potential :math:`\phi_{\text{Al}\text{Fe}}(r_{ij})`.
  * ``def ppfuncFeFe(r):`` - Fe-Fe pair-potential :math:`\phi_{\text{Fe}\text{Fe}}(r_{ij})`.

* **Embedding-Functions:**
  
  * ``def AlEmbedFunction(rho):`` - Al embedding function :math:`F_{\text{Al}}(\rho)`.
  * ``def FeEmbedFunction(rho):`` - Fe embedding function :math:`F_{\text{Fe}}(\rho)`.

* **Density-Functions:**
  
  * ``def AlAlDensityFunction(r):`` - Al density function :math:`\rho_{\text{Al}\text{Al}}(r_{ij})`.
  * ``def FeFeDensityFunction(r):`` - Fe density function :math:`\rho_{\text{Al}\text{Al}}(r_{ij})`.
  * ``def FeAlDensityFunction(r):`` - Al-Fe density function :math:`\rho_{\text{Al}\text{Fe}}(r_{ij})`.


.. note::

  The functional forms used within the Mendelev paper [#Mendelev2011]_ are somewhat long, and including their implementations here would detract from the readability of this example. However, they are included in the attached python file: :download:`eam_example3a.py`.


These functions are used within the ``main()`` function of the :download:`eam_example3a.py` file which is now shown:

.. literalinclude:: eam_example3a.py
  :language: python
  :pyobject: main

1. Breaking ``main()`` into its components, first a list of :class:`~atsim.potentials.Potential` objects is created, this is common with the other tabulation methods already discussed:

.. literalinclude:: eam_example3a.py
  :language: python
  :lines: 124-128

2. Next, the :class:`~atsim.potentials.EAMPotential` objects for Al and Fe are instantiated. This is where a Finnis-Sinclair model differs from the standard EAM seen earlier. Instead of a single density callable being passed to the :class:`~atsim.potentials.EAMPotential` constructor a dictionary of density functions is passed instead (see highlighted lines):
  
.. literalinclude:: eam_example3a.py
  :language: python
  :lines: 131-143
  :emphasize-lines: 4,5,10,11

3. The density function dictionary keys refer to the :math:`\beta` species in each :math:`\alpha\beta` pair. This means that:

  * for the Al :class:`~atsim.potentials.EAMPotential` instance:
    
    * :math:`\rho_{\text{Al}\text{Al}}` = ``AlAlDensityFunction()``,
    * :math:`\rho_{\text{Al}\text{Fe}}` = ``FeAlDensityFunction()``.
      
  * for the Fe :class:`~atsim.potentials.EAMPotential` instance:
    
    * :math:`\rho_{\text{Fe}\text{Al}}` = ``FeAlDensityFunction()``,
    * :math:`\rho_{\text{Fe}\text{Fe}}` = ``FeFeDensityFunction()``.


4. Finally, having defined the list of :class:`~atsim.potentials.EAMPotential` instances these are passed to the constructor of :class:`~atsim.potentials.eam_tabulation.SetFL_FS_EAMTabulation`\ , in this case writing the data to ``Mendelev_Al_Fe.eam.fs`` in the current directory:


.. literalinclude:: eam_example3a.py
  :language: python
  :lines: 152-156


Using the ``Mendelev_Al_Fe.eam.fs`` file within LAMMPS
""""""""""""""""""""""""""""""""""""""""""""""""""""""

For a binary system where Al and Fe have IDs of 1 and 2 the ``Mendelev_Al_Fe.eam.fs`` file is specified to LAMMPS as follows:

::
  
  pair_style eam/fs
  pair_coeff * * Mendelev_Al_Fe.eam.fs Al Fe



.. _eam_example_3b_obj:

Example 3b: Tabulate Al-Fe Finnis-Sinclair Potentials Using :class:`~atsim.potentials.eam_tabulation.TABEAM_FinnisSinclair_EAMTabulation` for DL_POLY
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Using exactly the same model definition as for :ref:`Example 3a <eam_example_3a_obj>`, the Al-Fe model can be re-tabulated for DL_POLY with minimal modification to the ``main()`` function. The modified version of the tabulation script can be found in :download:`eam_example3b.py`.  

The ``main()`` function is given below:

.. literalinclude:: eam_example3b.py
  :pyobject: main

Excluding the import statement at the top of the file, only two lines have been changed (highlighted). The first changes the filename to ``TABEAM`` whilst the second tells python to create an object of  :class:`~atsim.potentials.eam_tabulation.TABEAM_FinnisSinclair_EAMTabulation` instead of :class:`atsim.potentials.eam_tabulation.SetFL_FS_EAMTabulation`\ :

.. literalinclude:: eam_example3b.py
  :lines: 152-156
  :emphasize-lines: 1,4


That's it, nothing else has changed.


Using the ``TABEAM`` file with DL_POLY
""""""""""""""""""""""""""""""""""""""

Running :download:`eam_example3b.py` produces a file names ``TABEAM`` within the working directory. This should be placed in the same directory as the other DL_POLY input files (``CONTROL``, ``CONFIG`` and ``FIELD``). Then the following should be added to the end of the ``FIELD`` file:

::

  metal 3
  Al Al eeam
  Fe Fe eeam
  Al Fe eeam


.. note :: The Extended EAM (eeam) variant of the ``TABEAM`` file  generated here is only supported in DL_POLY versions >= 4.05.



.. rubric:: Footnotes

.. [#sutton1990] A.P. Sutton,  and  J. Chen,  "Long-range Finnis-Sinclair potentials", *Philos. Mag. Lett.* **61** (1990) 139 `doi:10.1080/09500839008206493 <https://dx.doi.org/10.1080/09500839008206493>`_\ .
.. [#Zhou2004] X. Zhou, R. Johnson and H. Wadley, "Misfit-energy-increasing dislocations in vapor-deposited CoFe/NiFe multilayers", *Phys. Rev. B.* **69** (2004) 144113.  
.. [#Mendelev2011] M.I. Mendelev, D.J. Srolovitz, G.J. Ackland, and S. Han, "Effect of Fe Segregation on the Migration of a Non-Symmetric Σ5 Tilt Grain Boundary in Al", *J. Mater. Res.* **20** (2011) 208.


.. _DL_POLY: http://www.stfc.ac.uk/cse/25526.aspx
.. _LAMMPS: http://lammps.sandia.gov

