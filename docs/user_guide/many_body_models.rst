.. _potable-many-body-models:

Many body models
================

Models that use the Embedded Atom Method (EAM) can be tabulated using :ref:`potable <potable-tool>`. Embedded atom models take the general form

.. math::
  :label: eq_standard_eam


  E_i = F_\alpha \left( \sum_{j \neq i} \rho_\beta(r_{ij}) \right) + \frac{1}{2} \sum_{j \neq i} \phi_{\alpha \beta} (r_{ij})


* Where:

  * :math:`\rho_\beta(r_{ij})` is the density function which gives the electron density for atom :math:`j` with species :math:`\beta` as a function of its separation from atom :math:`i`, :math:`r_{ij}`.
  * The electron density for atom :math:`i` is obtained by summing over the density (:math:`\rho_\beta (r_{ij}`) contributions due to its neighbours.
  * The embedding function :math:`F_\alpha(\rho)` is used to calculate the many-bodied energy contribution from this summed electron density.
  * The sum :math:`\frac{1}{2} \sum_{j \neq i} \phi_{\alpha \beta} (r_{ij})` gives the pair-potential contribution to atom :math:`i`'s energy. 
  * :math:`\phi_{\alpha \beta} (r_{ij})` are simply pair potentials that describe the energy between two atoms as a function of their separation.


In order to support the description of EAM when compared to  :ref:`pair-potential models <potable-pair-potential-models>`, additional sections in the input file are required. These are:

  * ``[EAM-Density]``: defines EAM density functions.
  * ``[EAM-Embed]``: describe the model's embedding functions.

As before, the pairwise section of the forcefield is specified in the :ref:`[Pair] section <potable-pair-section>` of the input file.

The EAM sections of the input are defined in much the same way as the :ref:`[Pair] section <potable-pair-section>`\ . Both  ``[EAM-Density]``\ , ``[EAM-Embed]`` allow the use of :ref:`multi-range potential definitions <multi-range-potentials>` and :ref:`potential modifiers <potential-modifiers>`\ . 

As for pair  only models, many-bodied force fields can specify :ref:`[Potential-Form] <potable-potential-form>` and  :ref:`[Table-Form] <potable-table-form>` sections can be provided if custom functional forms are required in ``[EAM-Density]``\ , ``[EAM-Embed]`` or ``[Pair]``\ .

``[EAM-Density]``
+++++++++++++++++

The density functions for embedded atom models are specified in this section. For the EAM form shown in equation :eq:`eq_standard_eam` entries in this section take the form::

    SPECIES : POTENTIAL_FORM PARAM_1 PARAM_2 ... PARAM_N


* Where:

  * ``SPECIES`` species for which density should be calculated
  * ``POTENTIAL_FORM PARAM_1 PARAM_2 ... PARAM_N`` potential form definition.

Density functions are tabulate in :math:`r_{ij}` space, therefore their extent and resolution are controlled by the the :ref:`dr <ref-potable-input-tabulation-dr>`\ , :ref:`nr <ref-potable-input-tabulation-nr>` and :ref:`cutoff <ref-potable-input-tabulation-cutoff>` fields in the :ref:`[Tabulation] <ref-potable-input-tabulation>` section in the same way as for ``[Pair]`` potentials.


Example
-------

A density function for silver may look something like this::

  [EAM-Density]
  Ag : as.exponential 4681.013008649 -6

This is taken from :ref:`potable_sutton_ag_example`\ . As you can determine from the parameters given there, and just for fun, this could have been defined using potential-modifiers as this, which makes the original parameters self evident::

  [EAM-Density]
  Ag : pow(
           product(as.constant 4.09, 
                   pow(
                       as.polynomial 0 1, 
                       as.constant -1)), 
           as.constant 6)

``[EAM-Embed]``
+++++++++++++++

Embedding functions are defined in this section. 

Entries have the following form::

    SPECIES : POTENTIAL_FORM PARAM_1 PARAM_2 ... PARAM_N

Where:

    * ``SPECIES`` is atomic type at which the surrounding electron density will be embedded using the specified potential form.
    * ``POTENTIAL_FORM PARAM_1 ...`` \ : embedding functions instantiate potential forms in the same way as in the :ref:`[Pair] <ref-potable-input-pair>` section.


.. note::

    Embedding functions are tabulated using rho values. The resolution and extent of functions in rho are defined by ``drho``\ , ``nrho`` and ``cutoff_rho`` in the :ref:`[Tabulation] section <ref-potable-input-tabulation>`\ .


Example
-------

This shows the embedding function used in the :ref:`potable_sutton_ag_example` for Ag::

  [EAM-Embed]
  Ag : product(as.constant 2.5415e-3, as.sqrt -144.41)


Note the use of the :ref:`product() <modifier-product>` modifier to apply a constant multiplication factor to the square root embedding function.


Standard EAM Examples
+++++++++++++++++++++

.. toctree::
  :hidden:

  potable_sutton_ag

More complete examples of EAM tabulation are listed in the following table:

================================ =======================================================================================
Example                          Description
================================ =======================================================================================
:ref:`potable_sutton_ag_example` An example of how to tabulate a single component EAM potential for Ag, to use in LAMMPS
================================ =======================================================================================


.. todo::

  Describe:

   * EAM Alloy
   * Finnis-Sinclair