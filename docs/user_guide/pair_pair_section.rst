.. _potable-pair-section:

******************
``[Pair]`` section
******************

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

.. seealso::

	See :ref:`list-of-potential-modifiers`\ .


.. _multi-range-potentials:

Multi-range potentials
++++++++++++++++++++++

The potential definition syntax used in the ``[Pair]`` section supports an extension which allows a series of potential-forms to be concatenated to each other, allowing each to act over a particular range of separations. These are defined as multi-range potentials. Concrete examples of where they are useful are provided in :ref:`multi-range-potential-examples` however the basic syntax defining multi-range potentials is introduced here. 

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


.. seealso::

	* Examples of multi-range potentials can be found here: :ref:`multi-range-potential-examples`