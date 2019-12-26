.. _potable-many-body-models:

Many body models
================

Models that use the Embedded Atom Method (EAM) can be tabulated using :ref:`potable <potable-tool>`. In comparison to :ref:`pair-potential models <potable-pair-potential-models>`, EAM models require additional sections in the input file. At a minimum it should contain:

  * ``[EAM-Density]``: defines EAM density functions.
  * ``[EAM-Embed]``: describe the model's embedding functions.
  * ``[Pair]``: the pair-potentials required by an EAM model.

.. _potable-eam-quick-start:

EAM Quick-Start
---------------

[1] S.M. Foiles, Application of the embedded-atom method to liquid transition metals, Phys. Rev. B. 32 (1985) 3409–3415.
10.1103/PhysRevB.32.3409

Density function:

.. math::

  \rho^a (r_{ij}) = (N-N_S) \rho_d^a(r_{ij}) + N_S \rho_S^a(r_{ij})

Where:
  * :math:`r_{ij}` is separation,
  * :math:`N` is the total number of outer electrons,
  * :math:`N_S` is a measure of the *s* electron content of the atomic density,
  * :math:`\rho_d^a` is the density of the outer *d* orbitals.

Embedding function:

The embedding energies are described by natural splines. 

.. code:: R

  cu_rho <- c(0.00000,0.01370,0.02740,0.05481,0.06303)
  cu_F <- c(0.00000,-2.92390,-4.29530,-2.82523, 0.00000)

  require(stats)

  f <- splinefun(cu_rho, cu_F, method="natural")
  ls(environment(f)$z)



Pair potential:

.. math::

  \Phi(r_{ij}) = Z^2(r_{ij})/r_{ij}


Where:

.. math::

  Z^2(r_{ij}) = a_1 (R_c - r_{ij})^3 + a_2 (R_c - r_{ij})^4

For :math:`r_{ij} < R_c` and zero when beyond :math:`R_c`.


Variants of the Embedded Atom Model
-----------------------------------

.. todo::

  Describe:

   * EAM Alloy
   * Finnis-Sinclair