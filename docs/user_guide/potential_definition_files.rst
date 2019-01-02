.. _using-potable:

***************
Using `potable`
***************

.. toctree:

	combining_potentials


This section of the user-guide will start by first explaining :ref:`pair-potential tabulation <potable-pair-potential-models>` before moving on to describe :ref:`many body potentials <potable-many-body-models>`, these contain additional terms that make their tabulation considerably more complicated.

.. contents::
	:local:



.. _potable-pair-potential-models:

Pair-potential models
=====================

The :ref:`Quick Start guide <quick_start>` should have already given you an idea of what a :ref:`potable <potable-tool>` potential definition looks like. We will now delve a little deeper and describe this input in more detail.

File Structure
--------------

Potable files for pair-potential models may contain the following sections::

	[Tabulation]
	...

	[Pair]
	...

	[Potential-Form]
	...

	[Species]
	...


Each section may contain a number of configuration options. These have the general form::

	ITEM : VALUE

or an equals sign may be used instead::

	ITEM = VALUE

Where ``ITEM`` identifies the configuration option's name and ``VALUE`` its value. 

Lines maybe commented out using ``#`` and line continuation is supported according to identation (see `here <https://docs.python.org/3/library/configparser.html#supported-ini-file-structure>`_ for more details). 


Each section of the file has a specific purpose and not all sections will be required in all cases:

* ``[Tabulation]`` describes how the file should be converted into a table file by the :ref:`potable <potable-tool` command. Contains information such as cutoff, output table format and cutoff. 
* ``[Pair]`` this is where pair interactions are defined by parametrising a potential-form.
* ``[Potential-Form]`` this section allows custom potential-forms to be defined. This may be required when you can't find an appropriate function from those :ref:`supplied <list-of-potential-forms>` with ``atsim.potentials``\ . However in many cases this won't be necessary and this section needn't appear in your model definition.
* ``[Species]`` this is used to provide meta-data about the species being tabulated. In most cases this section can be omitted (as very little species data is used during pair-tabulation). It is however, sometimes useful to include atomic charges etc, here so that the input file represents a complete description of a given potential model.


``[Tabulation]`` section
~~~~~~~~~~~~~~~~~~~~~~~~

This section defines the file format in which ``potable`` will write its output files through the :ref:`target <ref-potable-input-tabulation-target>` configuration item.

.. seealso::

	* Input format reference: :ref:`ref-potable-input-tabulation`\ .



The pair-tabulation :ref:`target <ref-potable-input-tabulation-target>` can be one of:

* ``DL_POLY`` (or ``DLPOLY``\ ): this creates output in the ``TABLE`` format accepted by the ``DL_POLY`` simulation code.
* ``LAMMPS``: creates output suitable for use with the LAMMPS ``pair_style table`` (see :ref:`using-table-in-lammps`\ ).
* ``GULP`` produces a table for the ``GULP`` code by defining a set of separation, energy pairs using the ``GULP`` ``spline`` directive.

Defining table cut-off
++++++++++++++++++++++

The extent of the table is defined in the ``[Tabulation]`` section using the :ref:`nr <ref-potable-input-tabulation-nr>`\ , :ref:`dr <ref-potable-input-tabulation-dr>` and :ref:`cutoff <ref-potable-input-tabulation-cutoff>` options:

* ``dr`` defines the row increment (step-size) between table rows.
* ``cutoff`` gives the maximum separation to be tabulated.
* ``nr`` determines the number of rows in the tabulation.


Any two of ``nr``\ , ``dr``\ and ``cutoff`` can be used to define the extent and resolution of the tabulation. As an example all three of the following ``[Tabulation]`` sections would produce a table with 1000 rows, each separated by 0.01:

#. ``cutoff`` and ``nr``\ ::

	[Tabulation]
	target: LAMMPS
	cutoff : 9.99
	nr : 1000


#. ``cutoff`` and ``dr``\ ::

	[Tabulation]
	target: LAMMPS
	cutoff : 9.99
	dr : 0.01


#. ``nr`` and ``dr``\ ::

	[Tabulation]
	target: LAMMPS
	nr : 1000
	dr : 0.01


.. _potable-pair-section:

``[Pair]`` section
~~~~~~~~~~~~~~~~~~
In the ``[Pair]`` section of the model definition, potential-forms are combined with parameters that tailor them to a given pair of species.

The *basic* form of an entry in this section is::

	SPECIES_A-SPECIES_B : POTENTIAL_FORM PARAM_1 PARAM_2 ... PARAM_N

* The label before the colon: ``SPECIES_A-SPECIES_B`` identifies the pair interaction being parametrised.
	+ This label consists of two species identifiers (``SPECIES_A`` and ``SPECIES_B``\ ) separated by a hyphen ``-``\ .
	+ The order in which the species labels are specified does not matter. That is, ``Au-Ag`` and ``Ag-Au`` would be equivalent.
	+ The ``[Pair]`` section can only contain one entry per unique species pair.
* The potential definition, specified after the colon, consists of the name of the potential-form (``POTENTIAL_FORM``\ ) followed by the numeric parameters it requires.

Pre-defined potential-forms
+++++++++++++++++++++++++++

A number of :ref:`pre-defined potential forms are provided  <list-of-potential-forms>`\ . These all have names pre-fixed by ``as.``

Each entry in the :ref:`list of potentials <list-of-potential-forms>` provides an entry called ``potable signature``\ . This shows the order in which parameters should be given to create a potential.

For the :ref:`Buckingham <potform-buck>` potential the ``potable signature`` is:

	| ``as.buck`` :math:`A` :math:`\rho` :math:`C`

which is associated with the formula:

.. math ::

	V(r_{ij}) = A \exp \left( - \frac{r_{ij}}{\rho} \right) - \frac{C}{r_{ij}^6}

This means that if we were defining a potential between Si and O that had :math:`A_{ij}` = 18003, :math:`\rho_{ij}` = 0.205 and :math:`C_{ij}` = 133.36 then the entry in the ``[Pair]`` section would be::

	[Pair]
	Si-O = as.buck 18003.0 0.205 133.36

Please refer to the ``potable signature`` when using the ``as.*`` potential-forms; specifying parameters in the wrong order will cause you problems.

It is also possible to define your own potential-forms in the ``[Potential-Form]`` section of ``potable`` file. These are parametrised here in the ``[Pair]`` section in the same way as the pre-defined ``as.*`` potential-forms. This usage is documented later here: :ref:`potable-potential-form`\ .


.. _potential-modifiers:

Potential modifiers
+++++++++++++++++++

If you followed the :ref:`quick_start` guide, you will have already seen a potential modifier. The ``[Pair]`` section from the :download:`basak.aspot <../quick_start/basak.aspot>` used in the :ref:`quick_start` is repeated here:

.. literalinclude:: ../quick_start/basak.aspot
	:lines: 6-10


You can see that the O-O and U-U pairs use the basic definition we have just seen. The U-O interaction however uses the modified form:


.. literalinclude:: ../quick_start/basak.aspot
	:lines: 9-10


Here ``sum()`` takes two basic pair-definitions (one for :ref:`as.buck <potform-buck>` and one for :ref:`as.morse <potform-morse>`\ ) and creates a pair-potential that is the sum of both. Here ``sum()`` is acting as a potential-modifier.

Potential-modifiers take the input or output of other potentials and produce outputs that have been altered in some way. A number of modifiers are provided with ``atsim.potentials`` and these are listed.

.. todo:: 
	
	Produce list of potential modifiers.


Multi-range potentials
++++++++++++++++++++++

The potential definition syntax used in the ``[Pair]`` section supports an extension which allows a series of potential-forms to be concatenated to each other, allowing each to act over a particular range of separations. These are defined as multi-range potentials. Concrete examples of where they are useful are provided in :ref:`combining-potentials` however the basic syntax defining multi-range potentials is introduced here. 

Suppose we want to define a potential acting between Mg and O using two potential-forms: ``pot_A`` and ``pot_B``.  The first is to be parametrised with values of 5.3 and 1.2 and ``pot_B`` with 9.6 and 2.4. Now say we want ``pot_A`` to act over the separations :math:`0 \geq r_{ij} \leq 3` and ``pot_B`` :math:`3 < r_{ij} \leq 8` and for the pair-potential to evaluate to zero when :math:`r_{ij} > 8`\ .

.. figure:: figures/multi_range.svg
	:name: fig_multi_range
	:width: 60%
	:align: center

	Illustration of multi-range potential definition described in the text.

The multi-range potential just described is summarised in :numref:`fig_multi_range`\ . This would be defined as follows in the ``potable`` input file::

	[Pair]
	Mg-O : >=0 pot_A 5.3 1.2 >3 pot_B 9.6 2.4 >8 as.zero


Notice that we used the :ref:`as.zero <potform-zero>` potential-form to provide a constant value of 0.0 when :math:`r_{ij} > 8` (equally ``as.constant 0.0`` could have been used). 

The syntax for a multi-range potential can be summarised as:

* A series of single potential definitions delimited by range markers.
* Range markers take the form:

    - ``>=R`` which indicates that the potential definition, following the marker, will be used at separations **greater than or equal** to the value specifid by ``R``
    - ``>R`` which means the same but acts only for separations **greater** than ``R``\ .


.. note::

	``potable`` defines all potentials to have the initial range of ``>0`` unless a range is explicitly defined. This is to avoid divide by zero errors when the potential is evaluated for :math:`r_{ij}` = 0. As this separation is unimportant to most physically relevant simulation.

	To include :math:`r_{ij}` = 0 in you tabulations simply make sure that your potential starts with ``>=0``\ ::

		>=0 POTENTIAL_DEFN


.. todo::

	Include see-also section for multi-range potentials.



``[Potential-Form]`` section
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This section of the input file is used for defining formulae that can be used as potential-forms and functions elsewhere in the model definition.

This allows for potential-forms that are not described in the standard file itself. Entries in this section have the general form::

	LABEL([ARG_1, ARG_2, ... , ARG_N]) : FORMULA


* ``LABEL`` is a unique identifier for the function (this can be used to refer to this function in the :ref:`[Pair] <potable-pair-section>`\ ).
* ``ARG ...`` defines the function signature by naming the arguments it takes.
* ``FORMULA`` mathematical expression defining the function.

Examples will now be given.

Example
+++++++

We will revisit the Basak [Basak2003]_ UO\ :sub:`2` model. In the :ref:`quick-start`\ guide you will have seen that the published potential parameters required considerable manipulation to make them compatible with the :ref:`potform-buck` and :ref:`potform-morse` potential-forms defined in ``atsim.potentials`` (see :ref:`basak-potential-model`\ ). Rather than transforming model parameters in this way, it may be easier to use the pair-potential equations and parameters directly as they appear in a paper. The ``[Potential-Form]`` section is the mechanism by which this may be achieved.

The Buckingham potential used in the Basak paper has the form:

.. math::

	V_\mathrm{Buck}(r_{ij}) & = f_0 b_{ij} \exp\left( \frac{a_{ij} - r_{ij}}{b_{ij}}\right) - \frac{C_{ij}}{r_{ij}^6} \\


We can write this as an entry in ``[Potential-Form]`` as::

	[Potential-Form]
	basak_buck(r,f0,a,b,c) = f0*b*exp((a-r)/b) - c/r^6


.. note::

	You may have noticed in this equation that we defined one of the terms using ``r^6`` (r to the sixth power), where python syntax would define this as ``r**6``\ . This is because, the formulae defined in this section are parsed using the `exprtk`_ library (via its python wrapper `cexprtk`_\ ). 

	To understand the functions, operators and syntax, supported for formulae please refer to the `exprtk`_ documentation. 


It is also possible to call the standard ``as.*`` :ref:`potential-forms <list-of-potential-forms>` in ``[Potential-Form]`` expressions. This is shown here to define a ``basak_morse`` function, where :ref:`as.morse <potform-morse>` will be used to provide a function that can be used directly with the parameters from the Basak paper (:numref:`table_basak_params`\ ). For reference, :ref:`as.morse <potform-morse> has the ``potable signature``\ :

	| ``as.morse`` :math:`\gamma` :math:`r_*` :math:`D`

Remembering that the :math:`D` parameter is given as :math:`f_0 \times D` using parameters from :numref:`table_basak_params` (see :ref:`basak-potential-model`\ ) we can now define our second potential function::

	[Potential-Form]
	basak_buck(r,f0,a,b,c) = f0*b*exp((a-r)/b) - c/r^6
	basak_morse(r, f0, d, gamma, r_star) = as.morse(r,gamma, r_star, f0*d)


.. note::

	It should also be noted that not-all the ``as.*`` potential-forms are available as functions within these formulae (for instance :ref:`as.buck4 <potform-buck4>` isn't). If you would like to check, please refer to the :ref:`list-of-potential-forms` and make sure that ``potential-function`` is listed as one of its ``Features``.


Now that we have both functions, we need to parametrise them for each interaction using values from :numref:`table_basak_params`\ . This is achieved in the normal way in the :ref:`[Pair] <potable-pair-section>`\ ::

	[Pair]
	O-O : basak_buck 0.042203 3.82 0.327022 3.948787
	U-U : basak_buck 0.042203 3.26 0.327022 0.0
	O-U : sum(
	        basak_buck  0.042203 3.54 0.327022 0.0,
	        basak_morse 0.042203 13.6765 1.65 2.369)


Notice that for the ``O-U`` interaction we continue to use the :ref:`sum() <modifier-sum>` to combine our Buckingham and Morse potentials (see :ref:`potential-modifiers`\ ).

The order in which parameters are specified in the ``[Pair]`` entries correspond to the arguments in the function signatures for ``basak_buck`` and ``basak_morse``, as is now shown:

.. figure:: figures/parameter_correspondence.*



.. note::

	By convention the ``as.*`` potentials take ``r`` (separation) as their first argument when used in formulae in the ``[Potential-Form]`` section.

	This represents a subtle to difference to when they appear in the ``[Pair]`` section and the argument list defined by the ``potable signature`` entries in :ref:`list-of-potential-forms`\ .

	For instance where ``as.buck`` could be parametrised as ``as.buck 1000.0 0.2 32.0`` in the ``[Pair]`` section it would be defined as ``as.buck(r, 1000.0, 0.2, 32.0)`` in a ``[Potential-Form]`` formula.


The model is now fully defined and gives the following potable input:

.. literalinclude:: example_files/basak_custom_potential_form_a.aspot


This input file can be downloaded as :download:`basak_custom_potential_form_a.aspot <example_files/basak_custom_potential_form_a.aspot>` and tabulated thus::

	potable basak_custom_potential_form_a.aspot Basak.lmptab

Section :ref:`using-table-in-lammps` describes how this table can then be used to perform a molecular dynamics simulation.


Alternative descriptions
++++++++++++++++++++++++

The potential-forms used in the previous example could have been defined in a number of different ways. Some of these are now shown to illustrate the flexibility of the ``potable`` system:

* :download:`basak_custom_potential_form_b.aspot <example_files/basak_custom_potential_form_b.aspot>`\ . In this example, a third potential-form ``basak_buckmorse`` is defined. This adds ``basak_buck()`` to ``basak_morse`` as an alternative to using the :ref:`sum() <modifier-sum>` potential modifier in the :ref:`[Pair] <potable-pair-section>`\ .

	.. literalinclude:: example_files/basak_custom_potential_form_b.aspot
	    :emphasize-lines: 15


* :download:`basak_custom_potential_form_c.aspot <example_files/basak_custom_potential_form_c.aspot>`. In this example the Morse potential is described directly rather than delegating to the``as.morse()`` function: 

	.. literalinclude:: example_files/basak_custom_potential_form_c.aspot
	    :emphasize-lines: 14

* :download:`basak_custom_potential_form_d.aspot <example_files/basak_custom_potential_form_d.aspot>`. This example shows that ``[Potential-Form]`` formulae can refer to each other.

    - In order to use the standard ``as.buck()`` potential-function, its :math:`A_{ij}` parameter must be calculated from the ``f0``\ , ``a`` and ``b`` Basak parameters (see :ref:`basak-potential-model`\ ).
    - Here an ``A_ij()`` formula is defined which is then invoked from inside the ``basak_buck()`` function. This sort of modularisation allows well structured and hence simpler expressions to be define.

	.. literalinclude:: example_files/basak_custom_potential_form_d.aspot
	    :emphasize-lines: 13,14




.. _potable-many-body-models:

Many body models
================


.. _exprtk: http://www.partow.net/programming/exprtk/index.html
.. _cexprtk: https://pypi.org/project/cexprtk/