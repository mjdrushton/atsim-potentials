.. _quick_start:

***********
Quick-Start
***********

The following provides a quick example of how a pair potential model can be defined and then tabulated for different simulation codes. For more advanced features such as EAM potentials, defining custom potential forms, splining or using the python interface, please see :ref:`user-guide`\ .

In this example we will define  Basak's [Basak2003]_  potential model for UO\ :sub:`2`\ . This has been chosen as it presents an issue that is often solved by using tabulated potentials. The U-O interaction combines Buckingham and Morse potential forms. Although most simulation codes natively provide these as analytical forms, some don't then allow them to be used in combination for the same pair-interaction. As a result, it becomes necessary to combine them externally and feed them into the code in tabulated form. This tutorial will show how this can be achieved using ``atsim.potentials``\ .

#. **Installation** (see :ref:`installation` for more detail)::

	pip install atsim.potentials

#. **Write input file**
	* In the Basak model [Basak2003]_  the O-O and U-U interactions are defined entirely using the Buckingham potential form. Whilst the O-U pair is the sum of a Buckingham and Morse potential (see :ref:`below <basak-potential-model>` for full details of the potential model). Both forms are provided by ``atsim.potentials`` and are specified as:

	.. math::
		:label: basak_short_components_standard_form

		V_\text{Buck}^\text{atsim}(r_{ij}) & = A_{ij} \exp\left( -\frac{r_{ij}}{\rho_{ij}} \right) - \frac{C_{ij}}{r_{ij}^6} \\
		V_\text{Morse}^\text{atsim}(r_{ij}) & = D_{ij} \left[ \exp \{ - 2 \gamma_{ij} (r_{ij} - r_{ij}^* \} - 2 \exp \{ - \gamma_{ij} (r_{ij} - r_{ij}^* \} \right] \\

	* The Buckingham potential is parametrised using values for :math:`A_ij`\ , :math:`\rho_{ij}` and :math:`C_{ij}`, specific to each species pair.
	* The Morse potential takes :math:`D_{ij}`\ , :math:`\gamma_{ij}`\ and :math:`r_{ij}^*`\ . 
	* Parameters for the Basak model are given in :numref:`table_basak_parameters_standard_form`\ .

	.. table:: Basak parameters for use with standard form of Buckingham and Morse potentials
		:name: table_basak_parameters_standard_form

		+---------------------------------------+-------------+------------+------------+
		|               Parameters              |     O-O     |    U-U     |    O-U     |
		+=======================================+=============+============+============+
		| :math:`A_{ij}`\ /eV                   | 1633.010243 | 294.640906 | 693.650934 |
		+---------------------------------------+-------------+------------+------------+
		| :math:`\rho_{ij}`\ /Å                 | 0.327022    | 0.327022   |   0.327022 |
		+---------------------------------------+-------------+------------+------------+
		| :math:`C_{ij}`\ /:math:`\text{eV}Å^6` | 3.948787    | 0.0        |        0.0 |
		+---------------------------------------+-------------+------------+------------+
		| :math:`D_{ij}`\ /eV                   | NA          | NA         |   0.577190 |
		+---------------------------------------+-------------+------------+------------+
		| :math:`\gamma_{ij}`\ /:math:`Å^{-1}`  | NA          | NA         |   1.650000 |
		+---------------------------------------+-------------+------------+------------+
		| :math:`r_{ij}^*`\ /Å                  | NA          | NA         |   2.369000 |
		+---------------------------------------+-------------+------------+------------+



	* The :ref:`potable <potable-tool>` tool is used to generate table files.
	* It accepts input in a straightforward format, reminiscent of ``.ini`` configuration files.
	* The potential parameters from :numref:`table_basak_parameters_standard_form` have been transferred into a file suitable for :ref:`potable <potable-tool>`\ : :download:`basak.aspot`\ :

		.. literalinclude:: basak.aspot

#. **Generate tabulated potential**
	* Download the :download:`basak.aspot` file and then generate a tabulation in ``DL_POLY`` format (see :ref:`below <specifying-other-tabulation-targets>` for information on other formats) by running this command::

		potable basak.aspot TABLE
	
	* This will generate a tabulation in the file named ``TABLE``\ .


And that's it! The rest of this page now gives details on what you've just done. Read on for more.

What are tabulated potentials?
==============================

Pair potentials are functions that relate potential energy to the separation of two interacting particles: :math:`U_{ij}(r_{ij})`\ . These are typically defined as equations, with a particular analytical form, that can be tailored to the chemistry of a given pair of species through a set of parameters (e.g. :math:`A_ij`\ , :math:`\rho_{ij}` and :math:`C_ij` in the case of the Buckingham form shown above).

Using these analytical potentials, an entire potential model, comprising of many interactions, can be described in a very compact form. Simulation codes come with large libraries of analytical potential forms, even so, it is impossible for all codes to support all potential forms. As a way of providing flexibility, and as an alternative to requiring users to edit and recompile simulation codes whenever they want to use an unusual form, most support table files.

Tabulated potential approximate a :math:`U_{ij}(r_{ij})` functions as a series of [separation, potential-energy] points that are stored in a file, readable by the code. For values not stored in the file, the simulation code performs interpolation to obtain energies and forces for in-between values.

:numref:`fig_tabulated_potential` shows the process followed in producing one of these tables. The desired analytical potential is defined and energy and forces (the first derivative of energy with respect to separation) are sampled at regular intervals. These samples are then written to a file in the format required by the simulation code.

The aim of ``atsim.potentials`` is to make this process simple, flexible and transferable across simulation codes. 

.. figure:: tabulated_potential.*
  :name: fig_tabulated_potential
  :width: 50%
  :align: center

  Tabulated potentials are obtained by sampling a mathematical formula at regular separations.


.. _basak-potential-model:

Potential model
===============

Before moving on to the tabulation procedure, it is useful to understand what it is we're trying to tabulate.


Basak's [Basak2003]_ model employs the following potential form:

.. math::
	:label: basak_eq1

	V(r_{ij}) = V_\mathrm{Coul}(r_{ij}) + V_\mathrm{Buck}(r_{ij}) + V_\mathrm{Morse}(r_{ij})


Where :math:`V(r_{ij})` is the potential energy of two atoms (:math:`i` and :math:`j`\ ) separated by :math:`r_{ij}`. The first term in eqn. :eq:`basak_eq1` defines the long range electrostatic interaction between the two atoms (with charges :math:`q_i` and :math:`q_j`):

.. math::
	:label: basak_coul

	V_\mathrm{Coul}(r_{ij}) = \frac{q_i q_j}{4\pi \epsilon_0 r_{ij}}

Where :math:`\epsilon_0` is the permittivity of free space. 

In a periodic system like UO\ :sub:`2`\ , it is necessary to employ various mathematical tricks to get the Coulomb sum to converge quickly and reliably (see for instance Ewald, cell-multipole or particle-mesh methods). As a result the Coulomb part of a potential model isn't normally included in a tabulation file and therefore it won't appear further in this example.

This leaves the short-range components of :math:`V(r_{ij})` from eqn. :eq:`basak_eq1` for us to bother about, namely :math:`V_\mathrm{Buck}(r_{ij})` and :math:`V_\mathrm{Morse}(r_{ij})`. In Basak's [Basak2003]_ paper these are defined as follows:

.. math::
	:label: basak_short_components

	V_\mathrm{Buck}(r_{ij}) & = f_0 b_{ij} \exp\left( \frac{a_{ij} - r_{ij}}{b_{ij}}\right) - \frac{C_{ij}}{r_{ij}^6} \\
	V_\mathrm{Morse}(r_{ij}) & = f_0 d_{ij} \left[ \exp \{ - 2 \gamma_{ij} (r_{ij} - r_{ij}^* \} - 2 \exp \{ - \gamma_{ij} (r_{ij} - r_{ij}^* \} \right] \\

The :math:`f_0`, :math:`a_{ij}`, :math:`b_{ij}`, :math:`C_{ij}`, :math:`D_{ij}`, :math:`\gamma_{ij}` and :math:`r_{ij}^*` specifying each pair interaction are given in the following table.

.. table:: Potential parameters for Basak model taken from original paper [Basak2003]_ (units have been converted for certain values).
  :name: table_basak_params

  +------------------------------------------+----------+----------+----------+
  | Parameters                               | O-O      | U-U      | O-U      |
  +==========================================+==========+==========+==========+
  | :math:`a_{ij}`\ /\ :math:`Å`             | 3.82     | 3.26     |     3.54 |
  +------------------------------------------+----------+----------+----------+
  | :math:`b_{ij}`\ /\ :math:`Å`             | 0.327022 | 0.327022 | 0.327022 |
  +------------------------------------------+----------+----------+----------+
  | :math:`C_{ij}`\ /\ :math:`\text{eV} Å^6` | 3.948787 | 0.0      |      0.0 |
  +------------------------------------------+----------+----------+----------+
  | :math:`d_{ij}`\ /\ :math:`r_{ij}^*`      | NA       | NA       |  13.6765 |
  +------------------------------------------+----------+----------+----------+
  | :math:`\gamma_{ij}`\ /\ :math:`r_{ij}^*` | NA       | NA       |     1.65 |
  +------------------------------------------+----------+----------+----------+
  | :math:`r_{ij}^*`                         | NA       | NA       |    2.369 |
  +------------------------------------------+----------+----------+----------+
  | :math:`f_0`\ /\ :math:`\text{eV} Å^{-1}` | 0.042203 | 0.042203 | 0.042203 |
  +------------------------------------------+----------+----------+----------+
  

The ``atsim.potentials`` package comes pre-loaded with a large number of potential-forms, so that you don't have to constantly redefine functional forms (see :ref:`list-of-potential-forms`). In this tutorial you will use the pre-defined versions of the Buckingham and Morse potentials.

If you compare the two definitions of the Buckingham and Morse potentials given in eqns. :eq:`basak_short_components_standard_form` and :eq:`basak_short_components` you will see there are some differences. The definitions of the Morse potential are almost identical, allowing :math:`\gamma_{ij}` and :math:`r_{ij}^*` to be used directly with the ``atsim`` Morse form. :math:`D_{ij}`, is however obtained as the product of :math:`f_0` and :math:`d_{ij}`

.. math::
	:label: basak_Dij

	D_ij = f_0 d_{ij}

The Buckingham potential used by Basak has more significant differences to the ``atsim`` standard. The :math:`C_{ij}` parameter can be used directly in both versions, however we need to manipulate the function from eqn. :eq:`basak_short_components` to show how :math:`A_{ij}` and :math:`\rho_{ij}` may be obtained from the original Basak parameter set (:numref:`table_basak_params`\ ).

Let's do this now by rearranging the Buckingham potential from eqn. :eq:`basak_short_components`:

.. math::
	:label: basak_buck_rearrange

	V_\mathrm{Buck}(r_{ij}) & = f_0 b_{ij}  \exp\left( \frac{a_{ij} - r_{ij}}{b_{ij}}\right) - \frac{C_{ij}}{r_{ij}^6} \\
													& = f_0 b_{ij}  \exp\left( \frac{a_{ij}}{b_{ij}} - \frac{r_{ij}}{b_{ij}} \right) - \frac{C_{ij}}{r_{ij}^6} \\
													& = f_0 b_{ij}  \frac{\exp\left( \frac{a_{ij}}{b_{ij}} \right)}{\exp \left( \frac{r_{ij}}{b_{ij}} \right)} - \frac{C_{ij}}{r_{ij}^6} \\
													& = f_0 b_{ij} \exp \left( \frac{ a_{ij} }{ b_{ij} } \right) \exp \left( -\frac{r_{ij}}{b_{ij}} \right) - \frac{C_{ij}}{r_{ij}^6} \\

By comparing coefficients between this equation, eqn. :eq:`basak_buck_rearrange` and the ``atsim`` form, eqn. :eq:`basak_short_components_standard_form`, it becomes obvious that the Basak parameters can be brought into the standard form using these relationships:

.. math::
	:label: basak_transformations

	A_{ij} & = f_0 b_{ij} \exp \left ( \frac{a_{ij}}{b_{ij}} \right ) \\
	\rho_{ij} & = b_{ij} \\

Using these relationships together with :math:`D_{ij}` from eqn. :eq:`basak_Dij` the table of potential parameters given above was obtained (:numref:`table_basak_parameters_standard_form`\ ).



Writing the potential definition
================================

Using the parameters from :numref:`table_basak_parameters_standard_form`\ , the model was described in the :download:`basak.aspot` file:

.. _basak_aspot:

.. literalinclude:: basak.aspot


This file contains two configuration blocks: 

* ``[Tabulation]``\ : this defines the table's output format

	+ ``target : DL_POLY`` means a DL_POLY TABLE file will be produced.
	+ ``cutoff : 6.5`` gives the maximum separation to include in the tabulation (:math:`r_{ij} = 6.5` Å).
	+ ``nr : 652`` the tabulation should contain 652 rows.

* ``[Pair]``\ : this section defines the O-O, U-U and O-U pair interactions:

	+ The basic form of each line is ``SPECIES_A-SPECIES_B = POTENTIAL_FORM PARAMETER_1 PARAMETER_2 ... PARAMETER_N``

	  - Where ``SPECIES_A`` and ``SPECIES_B`` define each pair of interacting species.
	  - ``POTENTIAL_FORM`` is a label identifying the functional form of the interaction. Here the ``as.buck`` and ``as.morse`` types are used. The ``as.`` prefix indicates these are standard forms provided by ``atsim.potentials`` (see :ref:`list-of-potential-forms` for a complete list).

The  :ref:`potform-buck` potential takes three parameters, separated by spaces::

	as.buck A ⍴ C


And the :ref:`potform-morse` parameters are::

	as.morse ɣ r* D

	
By comparing the :ref:`input file <basak_aspot>` and :numref:`table_basak_parameters_standard_form`, it should be apparent how the pair-interactions have been parametrised.

The O-U interaction makes use of a :ref:`potential modifier <potential_modifiers>` to combine the ``as.buck`` and ``as.morse`` forms. This is achieved using the :ref:`sum() <sum_modifier>` modifier which adds up all the contributions of the comma separated list of potentials defined inside the brackets.


Generating the tabulation
=========================

Save the :download:`basak.aspot` to your drive. Then generate the ``TABLE`` file by executing the :ref:`potable <potable-tool>` command::

	potable basak.aspot TABLE


This will interpret the potential model described above and write it in ``DL_POLY`` format to a file named ``TABLE``.

For an example of how this ``TABLE`` file can be used in a ``DL_POLY`` simulation see :ref:`using-table-in-dlpoly`.

.. toctree::
  :hidden:

  ../basak_tabulate_dlpoly/using_table_in_dlpoly


.. _specifying-other-tabulation-targets:

Specifying other tabulation targets
===================================

Once the potential model has been defined, creating tabulations for different codes in different formats is simple. 

LAMMPS
------

To target LAMMPS it is just a case of changing the ``target`` option in the ``[Tabulation]`` section of the :download:`basak.aspot` file to ``LAMMPS`` e.g. ::

	[Tabulation]
	target : LAMMPS
	cutoff : 6.5
	nr : 652


A ``LAMMPS`` table file named ``Basak.lmptab`` would then be generated by re-running :ref:`potable <potable-tool>`::

	potable basak.aspot Basak.lmptab


An example of how to use this file in LAMMPS is given here: :ref:`using-table-in-lammps`.

GULP
----

In the previous section a new tabulation target was specified by editing the :download:`basak.aspot` file. For temporary changes, the :ref:`potable <potable-tool>` allows configuration options to be overridden from the command line. Let's do this now to make a table file suitable for GULP. This is achieved with the :ref:`--override-item <cmdoption-override-item>` option, like this::

	potable --override-item Tabulation:target=GULP basak.aspot potentials.lib


The ``potentials.lib`` file can be used with the following GULP input file to run an energy minimisation: :download:`basak.gin`


.. toctree::
  :hidden:

  ../basak_tabulate_lammps/using_table_in_lammps


 
.. _DL_POLY: http://www.stfc.ac.uk/cse/25526.aspx
