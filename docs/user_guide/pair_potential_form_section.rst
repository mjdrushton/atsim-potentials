.. _potable-potential-form:

****************************
``[Potential-Form]`` section
****************************

This section of the input file is used for defining formulae that can be used as potential-forms and functions elsewhere in the model definition.

This allows for potential-forms that are not described in the standard file itself. Entries in this section have the general form::

	LABEL([ARG_1, ARG_2, ... , ARG_N]) : FORMULA


* ``LABEL`` is a unique identifier for the function (this can be used to refer to this function in the :ref:`[Pair] <potable-pair-section>`\ ).
* ``ARG ...`` defines the function signature by naming the arguments it takes.
* ``FORMULA`` mathematical expression defining the function.

Examples will now be given.

.. _potable-potential-form-basak-example:

Example: using custom-potential forms to define Basak potential
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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


Now that we have both functions, we need to parametrise them for each interaction using values from :numref:`table_basak_params`\ . This is achieved in the normal way in the :ref:`[Pair] <potable-pair-section>` section::

	[Pair]
	O-O : basak_buck 0.042203 3.82 0.327022 3.948787
	U-U : basak_buck 0.042203 3.26 0.327022 0.0
	O-U : sum(
	        basak_buck  0.042203 3.54 0.327022 0.0,
	        basak_morse 0.042203 13.6765 1.65 2.369)


Notice that for the ``O-U`` interaction we continue to use the :ref:`sum() <modifier-sum>` potential-modifier to combine our Buckingham and Morse potentials (see :ref:`potential-modifiers`\ ).

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



.. _exprtk: http://www.partow.net/programming/exprtk/index.html
.. _cexprtk: https://pypi.org/project/cexprtk/