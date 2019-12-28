.. _potable-many-body-models:

Many body models
================

Models that use the Embedded Atom Method (EAM) can be tabulated using :ref:`potable <potable-tool>`. In comparison to :ref:`pair-potential models <potable-pair-potential-models>`, EAM models require additional sections in the input file. At a minimum it should contain:

  * ``[EAM-Density]``: defines EAM density functions.
  * ``[EAM-Embed]``: describe the model's embedding functions.
  * ``[Pair]``: the pair-potentials required by an EAM model.

As for pair  only models, many-bodied force fields can specify :ref:`[Potential-Form] <potable-potential-form>` and  :ref:`[Table-Form] <potable-table-form>` sections can be provided if custom functional forms are required.

.. seealso::

  * **Reference guide:** :ref:`ref-potable-eam-density`




Variants of the Embedded Atom Model
-----------------------------------

.. todo::

  Describe:

   * EAM Alloy
   * Finnis-Sinclair