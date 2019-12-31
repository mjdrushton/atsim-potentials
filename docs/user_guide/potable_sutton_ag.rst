

.. _potable_sutton_ag_example:

Sutton Ag EAM Example
=====================

This provides an example of using ``potable`` to tabulate the Ag model given by Sutton and Chen in [#sutton1990]_\ .

Potential Model
+++++++++++++++

.. math::
    :label: eq_sutton_model

    E_i = \epsilon \left[ \frac{1}{2} \underset{i \neq j}{\sum \sum} \phi_{\alpha\beta}(r_{ij}) - c \sum_i \sqrt{\rho_i} \right]

Where: 

    * :math:`E_T` is the energy of ths system.
    * :math:`r_{ij}` is the separation between atoms :math:`i` and :math:`j`\ .
    * :math:`c` and :math:`\epsilon` are adjustable parameters specific to interacting species.
    * Inside the square brackets the first term :math:`V(r_{ij})` are the pair potentials.
    * The second is the many body term: :math:`c \sum_i \sqrt{\rho_i}`\ . Where :math:`\rho_i` is the electron density.

Pair potential form:
--------------------

.. math::

    \phi_{\alpha \beta}(r_{ij}) = (a/r_{ij})^n

Where: 

    * :math:`a` and :math:`n` are potential parameters.

This must be multiplied by the :math:`\epsilon` term from equation :eq:`eq_sutton_model` above:

.. math::

    \phi_{\alpha \beta}(r_{ij}) = \epsilon (a/r_{ij})^n


To make things easier later on, this will be re-expressed as:

.. math::

    \phi_{\alpha \beta}(r_{ij}) = \epsilon a^n r_{ij}^n


This will allow this functional form to be written using the provided :ref:`as.exponential <potform-exponential>` potential-form.

Many body terms
---------------

Density function:
^^^^^^^^^^^^^^^^^

The density function is:

.. math::

    \rho_i = \left( \frac{a}{r_{ij}} \right)^m

Again to allow the use of the :ref:`as.exponential <potform-exponential>` potential-form this will be re-written as:

.. math::

    \rho_i = a^m r_{ij}^{-m}


Embedding function:
^^^^^^^^^^^^^^^^^^^

Examining the many-body term from :eq:`eq_sutton_model` it can be seen that the embedding function is:

.. math::

    c \sqrt{\rho_i}

Taking the the :math:`\epsilon` term from outside the square brackets and pre-multiplying the expression this becomes:

.. math::

    \epsilon c \sqrt{\rho_i}


Potential parameters
--------------------

The potential parameters for Ag are:

.. table:: Potential parameters for Ag
    :name: tab_potable_sutton_ag_potential_parameters

    =================  =======================
    Parameter          Value
    =================  =======================
    :math:`m`          6
    :math:`n`          12
    :math:`\epsilon`   2.5415×10\ :sup:`-3` eV
    :math:`a`          4.09
    :math:`c`          144.41
    =================  =======================


Potable input
+++++++++++++

The ``potable`` input for this model can be downloaded as :download:`Ag_sutton.aspot <example_files/Ag_sutton.aspot>` and will now be described:

.. literalinclude:: example_files/Ag_sutton.aspot
    :linenos:

.. rubric:: Notes:

* **lines 1-8 [Tabulation]:** 
    
    * **lines 4,5:** gives the resolution and extent of the function in ``[EAM-Embed]``\ .
    * **lines 7,8:** defines resolution and extent of the tables generated for the ``[Pair]`` and ``[EAM-Density]`` functions.

* **lines 10 and 11 [EAM-Embed]:**
    
    * Defines the embedding function.
    * Note the use of the ``product()`` :ref:`potential modifier <potential-modifiers>` to multiply the square root embedding function by the value of :math:`\epsilon`\ .

* **lines 13 and 14 [EAM-Density]:**

    * Describes the density function.
    * The value of 4681.013008649 is obtained as :math:`a^m = 4.09^6`\ .

* **lines 16 and 17 [Pair]:**

    * Defines the pair potential component of the model.
    * As above, the ``product()`` :ref:`potential modifier <potential-modifiers>` has been used to multiply the function by :math:`\epsilon`\ .
    * Here the first parameter to the ``as.exponential`` form is :math:`a^n = 4.09^{12}` = 21911882.78.

Making and testing the tabulation
+++++++++++++++++++++++++++++++++

To tabulate the potential :download:`download the aspot file <example_files/Ag_sutton.aspot>` and run it through ``potable``\ ::

    potable Ag_sutton.aspot Ag_sutton.eam.alloy



A LAMMPS input file is provided to allow you to test the ``Ag_sutton.eam.alloy`` file produced by potable. This input file can be downloaded here: :download:`Ag_sutton_fcc.lmpin <example_files/Ag_sutton_fcc.lmpin>` and will energy minimize the structure and then perform an NPT MD equilibration at T=300K. Frames will be dumped every 1000 timesteps (1ps) and dumped to a LAMMPS dump file named ``dump.atom``\ (this is suitable for visualisation in `Ovito <https://ovito.org>`_\ ).

In terms of the table file the important part of the LAMMPS input is::

    pair_style eam/alloy
    pair_coeff * * Ag_sutton.eam.alloy Ag


This tells LAMMPS to accept a ``setfl`` formatted file (``pair_style eam/alloy``\ ). The ``Ag`` at the end of the ``pair_coeff`` line says that LAMMPS should associate atom type 1 with the ``Al`` species label in the table file ``Ag_sutton_eam.alloy``\ .


Placing both the LAMMPS and table file in the same directory run LAMMPS as follows::

    mpirun lammps -in Ag_sutton_fcc.lmpin -log Ag_sutton_fcc.lmpout


.. rubric:: Footnotes:

.. [#sutton1990] A.P. Sutton,  and  J. Chen,  "Long-range Finnis-Sinclair potentials", *Philos. Mag. Lett.* **61** (1990) 139 `doi:10.1080/09500839008206493 <https://dx.doi.org/10.1080/09500839008206493>`_\ .