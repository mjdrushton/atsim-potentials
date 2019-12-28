.. _potable-many-body-models:

Many body models
================

Models that use the Embedded Atom Method (EAM) can be tabulated using :ref:`potable <potable-tool>`. In comparison to :ref:`pair-potential models <potable-pair-potential-models>`, EAM models require additional sections in the input file. At a minimum it should contain:

  * ``[EAM-Density]``: defines EAM density functions.
  * ``[EAM-Embed]``: describe the model's embedding functions.
  * ``[Pair]``: the pair-potentials required by an EAM model.

As for pair  only models, many-bodied force fields can specify :ref:`[Potential-Form] <potable-potential-form>` and  :ref:`[Table-Form] <potable-table-form>` sections can be provided if custom functional forms are required.


Example: Voter-Chen Aluminium and Nickel Embedded Atom Model
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The following example describes how to tabulate the Al and Ni potential models developed by Voter and Chen and described in the following paper:

  * A.F. Voter and S.P Chen,  "Accurate Interatomic Potentials for Ni, Al and Ni\ :sub:`3`\ Al",  *MRS Proceedings*\ ,  `doi:10.1557/proc-82-175 <https://dx.doi.org/doi:10.1557/proc-82-175>`_ (1986) 82.

This example will show how to generate a ``setfl`` embedded atom model tabulation suitable for use with the LAMMPS `eam/alloy <https://lammps.sandia.gov/doc/pair_eam.html#pair-style-eam-alloy-command>`_ pair_style.


.. _potable-many-body-models-voter-chen-potential-form:

Potential Form
--------------

The different components of the potential model will now be described. 

**Pair:** The paper uses the Morse pair potential form. Whilst we could use the predefined form provided with ``atsim.potentials``, for consistency with the paper we use the following:

.. math::

  \phi (r)  = D_M \left\{ 1- \exp \left[  - \alpha_M (r-R_M) \right] \right\}^2 - D_M

Where :math:`D_M`\ , :math:`\alpha_M`\  and :math:`R_M` are variable parameters.

**Density function:** The following density function was used which is based on a hydrogenic 4s orbital. Here :math:`\beta` is an adjustable parameter:

.. math::

  \rho (r) = r^6 \left[ \exp(-\beta r) + 2^9 \exp(-2 \beta r) \right]

**Embedding function:** Voter and Chen's paper used a method to determine the embedding function which is not well suited to an analytical definition. Therefore it was digitised from the paper using `Web Plot Digitizer <https://automeris.io/WebPlotDigitizer/>`_ and will be specified as a :ref:`[Table-Form] <potable-table-form>` in the ``potable`` input below (see :numref:`fig_voter_chen_embedding_functions`\ ).

.. figure:: figures/voter_chen_embedding_functions.*
  :name: fig_voter_chen_embedding_functions

  The aluminium and nickel embedding functions digitised from the Voter and Chen paper.


Writing the ``potable`` files
-----------------------------

As both use the same density and pair-potential forms (see :ref:`above <potable-many-body-models-voter-chen-potential-form>`\ ) let us define these first::

  [Potential-Form]
  morse(r, D_M, R_M, alpha_M) = D_M * (1 - exp(-alpha_M * (r-R_M)))^2 - D_M
  density(r, beta) = r^6 * ( exp(-beta * r) + 2^9 * exp(-2*beta*r))

:numref:`table_voter_chen_potential_parameters` provides potential parameters that will now be used to write two ``potable`` files, one for each metal. 

.. table:: Voter-Chen potential parameters for Ni and Al. Pair parameters are :math:`D_M` :math:`R_M`\ , :math:`\alpha_M` and :math:`r_\text{cut}`\ . The  :math:`\beta` parameter is relevant to the EAM density function.
  :name: table_voter_chen_potential_parameters

  =================================  ====== ======
  Potential parameter                Ni     Al    
  =================================  ====== ======
  :math:`D_M` (eV)                   1.5335 3.7760
  :math:`R_M` (Å)                    2.2053 2.1176
  :math:`\alpha_M` (Å\ :sup:`-1`\ )  1.7728 1.4859  
  :math:`r_\text{cut}` (Å)           4.7895 5.5550
  :math:`\beta` (Å\ :sup:`-1`\ )     3.6408 3.3232
  =================================  ====== ======


Aluminium
~~~~~~~~~

The Voter-Chen tabulation for aluminium can be downloaded as: :download:`Al_voter_chen.aspot <example_files/Al_voter_chen.aspot>`\ .

The :download:`Al_voter_chen.aspot <example_files/Al_voter_chen.aspot>` file can be used to generate a ``setfl`` file (``Al_voter_chen.eam.alloy``\ ) by running ``potable``\ ::

  potable Al_voter_chen.aspot Al_voter_chen.eam.alloy


A LAMMPS input file is provided to use with the ``Al_voter_chen.eam.alloy`` file created by this command. This can be downloaded here: :download:`Al_voter_chen_equil.lmpin <example_files/Al_voter_chen_equil.lmpin>`\ . This defines a 4×4×4 Al FCC super-cell then performs an NPT equilibration at T=300K. With the LAMMPS input and table file in the same directory this could be run as follows::

  mpirun lammps -in Al_voter_chen_equil.lmpin

The following lines in the LAMMPS input make use of the setfl table file::

  pair_style eam/alloy
  pair_coeff * * Al_voter_chen.eam.alloy Al

The trailing ``Al`` states that atom type 1 is identiifed by the label ``Al`` in the table file.


The :download:`potable input file <example_files/Al_voter_chen_equil.lmpin>` has the following contents:

.. literalinclude:: example_files/Al_voter_chen.aspot
  :linenos:


**Notes:**

  * **line 2:** ``target : setfl`` showing this is an EAM tabulation.
  * **lines 3 and 4:**  additional step and cutoff values are defined for embedding function which works in ``rho`` space.
  * **lines 8 and 9:** parameterises the ``morse`` potential using the definition on line 18. It specifies the Al-Al pair interaction's :math:`D_M`\ , :math:`R_M` and :math:`\alpha_M` parameters using the values given in :numref:`table_voter_chen_potential_parameters`\ . Also note that this is a :ref:`multi-range potential <multi-range-potentials>` that becomes zero above the value of :math:`r_\text{cut}` given for Al (5.5550 Å)::

      [Pair]
      Al-Al : morse 3.7760 2.1176 1.4859 >=5.5550 as.zero


  * **lines 11-15:** provides the details of the EAM density and embedding function:
  
    * **[EAM-Density] lines 11-12:** this creates the Al density function by parameterising the functional form from line 19 with a :math:`\beta` value of 3.3232::
    
        [EAM-Density]
        Al : density 3.3232


    * **[EAM-Embed] lines 11-12:** here the :ref:`[Table-Form] <potable-table-form>` described between lines 21 and 51 (named ``Al_embed``\ ) is specified as the Al embedding function. Notice that this potential-form does not take any parameters::

        [EAM-Embed]
        Al : Al_embed



Nickel
~~~~~~

The ``potable`` file for Ni is very similar as can be seen below. Similarly a LAMMPS equilibration file is also provided here: :download:`Ni_voter_chen_equil.lmpin <example_files/Ni_voter_chen_equil.lmpin>`\ .

The potable input :download:`(Ni_voter_chen_equil.aspot) <example_files/Ni_voter_chen_equil.aspot>` is now shown:

.. literalinclude:: example_files/Ni_voter_chen_equil.aspot









Variants of the Embedded Atom Model
-----------------------------------

.. todo::

  Describe:

   * EAM Alloy
   * Finnis-Sinclair