.. _combining-potentials:

********************
Combining potentials
********************

It is often useful to describe the different regions of a pair-potential using different functional forms. This section of the user-guide will show ways in which different potentials can be combined.

First, examples will be given demonstrating the use :ref:`multi-range potentials <multi-range-potentials>` described earlier.

Methods for linking potential-forms using splines will then be considered.


.. _multi-range-potential-examples:

Multi-range potential examples
==============================

This section shows how potential-forms can be combined using the multi-range syntax introduced here: :ref:`multi-range-potentials`\ .

Example: truncating a potential
-------------------------------

The starting structure for a simulation may contain overlapping atoms (for example if generated from random coordinates). In these cases it is useful to move atoms apart before starting the simulation proper.

In principle this should be possible by using an energy minimisation run, however some potential forms are badly behaved at very small separations (for instance the Buckingham catastrophe) leading to unexpected results.

One way of removing atom overlap is to temporarily substitute the potentials with some that are guaranteed to work, even when two atoms are completely on top of each other. An example of a potential-form that could be used for this purpose is now described (this is available natively in LAMMPS as `pair_style soft`_\ ):

.. math::

	V_{ij}(r_{ij}) = A \left[ 1+ \cos \left( \frac{\pi r_{ij}}{r_c} \right) \right]


Let's start by plotting this function with sensible parameter values of :math:`A` = 10 and :math:`r_c` = 1.6:


.. figure:: figures/soft_plot.*
	:width: 60%
	:align: center


To make this useful, the potential is truncated at :math:`r_c` (i.e. :math:`r_{ij} < r_c`\ ). This means that the potential is repulsive below :math:`r_c` and doesn't act above:

.. figure:: figures/soft_plot_truncated.*
	:width: 60%
	:align: center


Now, having been suitably truncated, the potential will gently push atoms apart in an MD or energy minimisation run. This can be implemented for :ref:`potable <potable-tool>` as follows; first define the cosine function:

.. literalinclude:: example_files/soft_a.aspot
	:lines: 8,9


The truncation at :math:`r_c` can then be applied using a multi-range ``[Pair]`` deinition. 

.. literalinclude:: example_files/soft_a.aspot
	:lines: 6-8


Here the Si-O interaction have :math:`A` = 10.0 and an :math:`r_c` value of 1.6Å (the position of the Si-O peak in a silica glass radial distribution function) after which the :ref:`potform-zero` potential-form takes over. Similarly, the O-O pair has a weaker repulsion (:math:`A`=5.0 \ ) and a cutoff-radius, :math:`r_c`\ , of 2.4Å (again this is about the position of the O-O peak in the RDF for silica glass).

In this way the cosine function has been truncated. Specifying an Si-Si interaction will be left as an exercise for the reader.

The complete input file for this example can be downloaded here :download:`soft_a.aspot <example_files/soft_a.aspot>` and is:

.. literalinclude:: example_files/soft_a.aspot



Example: truncating using ``if()`` statement
--------------------------------------------

This example shows an alternative to using the multi-range syntax to implement the soft cosine potential described above. The `exprtk`_ package used to evaluate expressions in the :ref:`[Potential-Form] <potable-potential-form>` section, support an ``if ... then ... else`` construct through its ``if()`` function. This has the form::

	if (CONDITION, FORMULA_IF_TRUE, FORMULA_IF_FALSE)

This means the :math:`r > r_c` condition could have been included in :ref:`[Potential-Form] <potable-potential-form>` like this:

.. literalinclude:: example_files/soft_if.aspot
	:lines: 10-12


This simplifies the definition of the ``[Pair]`` section by avoiding the repetition of the :math:`r_c` parameter when defining the multi-range version of the potential:

.. literalinclude:: example_files/soft_if.aspot
	:lines: 6-8


The :ref:`potable <potable-tool>` input using this version of the potential-form can be downloaded here: :download:`soft_if.aspot <example_files/soft_if.aspot>`\ .


.. _potable-published-spline-coefficients:

Example: parametrising a model using published spline coefficients
------------------------------------------------------------------

The UO\ :sub:`2` potential model by Morelon describes U-O and O-O interactions [Morelon2003]_\ . The former employs the :ref:`potform-bornmayer` potential-form while the latter uses the combination of a 3rd and 5th order polynomial spline to link an :math:`A \exp\left( - \frac{r_{ij}}{\rho}\right)` term to a :math:`-\frac{C}{r_{ij}^6}` term. The spline-coefficients for the Morelon model are available in a paper by Potashnikov et al. [Potashnikov2011]_\ . They give the O-O interaction as:

.. math::
	:label: eq_morelon-oo

	V(r_{ij}) = 
	\begin{cases}
		11272.6 \exp\left(\frac{0.1363}{r_{ij}}\right) &: r_{ij} < 1.2 \\
		479.955
		-1372.53 r_{ij}
		+1562.22 r_{ij}^2
		-881.969 r_{ij}^3
		+246.435 r_{ij}^4
		-27.2447 r_{ij}^5 &: 1.2 \leq r_{ij} < 2.1 \\
		 42.8917
		-55.4965 r_{ij}
		+23.0774 r_{ij}^2
		-31.13140 r_{ij}^3 &: 2.1 \leq r_{ij} < 2.6 \\
		-\frac{134}{r_{ij}^6} &: r_{ij} \geq 2.6 \\
	\end{cases}



The following figure plots the O-O interaction to show how its constituent parts relate to each other.


.. figure:: figures/morelon.svg
	:width: 60%
	:align: center


	The O-O interaction of the Morelon model is made up of several potential-forms acting over different ranges. The upper plot shows the un-trimmed versions of the individual functions. The lower plot shows how these combine to form the overall interaction. 

	**Note:** in a simulation the O-O interaction would also include a repulsive electrostatic term which isn't plotted here.


Using this information we can write the Morelon model as :ref:`potable <potable-tool>` input (and can be downloaded as :download:`morelon.aspot <example_files/morelon.aspot>`\ ) :

.. literalinclude:: example_files/morelon.aspot


There are a few features here that are worth noting:

* The four distinct regions are defined a multi-range potentials.
* The third and fifth order polynomials are both defined using the :ref:`as.polynomial <potform-polynomial>` form:

    - The order of the polynomial is implicitly defined by the number of parameters.
* The :math:`-\frac{C}{r_{ij}^6}` term is implemented using the :ref:`potform-buck` form. In order to only use the attractive part of the potential an ``A`` parameter of 0.0 is used as its first argument and a non-zero ``rho`` is then specified.

.. _aspot-splining:

Splining
========

Splines are curves that are used to smoothly connect different functional forms across different ranges of a potential-form. Typically splining curves are polynomials. Care must be taken when determining the coefficients of the spline function as it is important that potentials are smooth and do not have steps in their derivatives. For this reason spline parameters are determined so that at its end-points, it gives the same potential-energy and derivatives (1st and 2nd) as the functions it joins. 

In the previous example (see :ref:`potable-published-spline-coefficients`\ ), spline coefficients were provided but their calculation can be quite involved. As a result ``potable`` provides services to make this easier. 

Splining performed in the ``[Pair]`` section of the input through the ``spline()`` :ref:`potential-modifier <potential-modifiers>`\ . 


``spline()`` takes a single argument which has the form of a :ref:`multi-range potential definition <multi-range-potentials>`\ . It has three distinct parts, representing the three ranges of the splined potential:

#. The region described by the starting potential.
#. The interpolation region given by the spline.
#. The final part defined by the end potential.

In the description that follows, the separation at which the starting potential ends will be referred to as the **detachment-point**\ . Similarly, the end of the splined region is the **attachment-point**\ .

Bearing in mind its similarity to the :ref:`multi-range <multi-range-potentials>` syntax ``spline()`` definitions have this basic form::

	spline(POT_A_DEFN >R_DETACH SPLINE_DEFN >R_ATTACH POT_B_DEFN)

Where:

* ``POT_A_DEFN`` - this is the potential definition for the starting-potential. This is typically a potential-form label followed by potential parameters (see :ref:`potable-pair-section`\ ).
* ``R_DETACH`` -  detachment-point (**note:** as this is a multi-range potential definition the inclusive form ``>=R_DETACH`` is also valid).
* ``SPLINE_DEFN`` - the type of spline to be used along with any parameters it takes. This will be described properly in the next paragraph.
* ``R_ATTACH`` -  attachment point (again ``>=R_ATTACH`` is also accepted).
* ``POT_B_DEFN`` - the potential definition of the end-potential.

``SPLINE_DEFN`` has the same format as a potential-definition: an identifying label followed by a list of parameters (**note:** some spline-types do not take any parameters)::

	``SPLINE_LABEL PARAM_1 PARAM_2 ... PARAM_N``

Currently ``SPLINE_LABEL`` can be :ref:`exp_spline <spline-exp_spline>` or :ref:`buck4_spline <spline-buck4>`\ . These spline-types will now be described.


.. _spline-exp_spline:

Exponential Spline ``exp_spline``
---------------------------------

The spline used when ``exp_spline`` is specified is given as:

.. math::
	:name: eq_exp_spline

    
    V(r_{ij}) = \exp \left( B_0 + B_1 r_{ij} + B_2 r_{ij}^2 + B_3 r_{ij}^3 + B_4 r_{ij}^4 + B_5 r_{ij}^5 \right)


Where :math:`B_{0..5}` are the spline coefficients calculated automatically during splining. 

The exponential spline does not take any parameters.

The functional form given in eqn. :eq:`eq_exp_spline` is also available as a potential-form. Although the user is then responsible for calculating their own spline-coefficients. See: :ref:`potform-exp_spline`\ .

.. _spline-exp_spline-example:

Example: splining to the ``zbl`` potential form using ``exp_spline``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example will show how to use ``spline()`` with the ``exp_spline`` type to add the repulsive :ref:`potform-zbl` to a :ref:`potform-buck` potential that is unphysical at small separations.

For certain parameterisations, popular potential forms can exhibit unphysical behaviour for some interatomic separations. A popular model for the description of silicate and phosphate systems is that due to van Beest, Kramer and van Santen (the BKS potential set) [VanBeest1990]_\ . In the current example, the Si-O interaction from this model will be considered. This uses the Buckingham potential form with the following parameters:
    
    * A = 18003.7572 eV
    * :math:`\rho` = 0.205204 Å
    * C = 133.5381 eV :math:`Å^6`
    * Charges:
        -   Si = 2.4 *e*
        -   O  = -1.2 *e*

The plot in :numref:`fig_exp_spline_problem` shows the combined coulomb and short-range contributions for this interaction plotted as a function of separation. The large C term necessary to describe the equilibrium properties of silicates means that as :math:`r_{ij}` gets smaller, the :math:`\frac{C}{r_{ij}^6}` overwhelms the repulsive Born-Mayer component of the Buckingham potential meaning that it turns over. This creates only relatively shallow minimum around the equilibrium Si-O separation. Within simulations containing high velocities (e.g. high temperatures or collision cascades) atoms could easily enter the very negative, attractive portion of the potential at low :math:`r_{ij}` - effectively allowing atoms to collapse onto each other. In order to overcome this deficiency a ZBL potential will be splined onto the Si-O interaction within this example.

.. figure:: figures/exp_spline_problem.*
	:name: fig_exp_spline_problem
	:width: 60%
	:align: center

	Plot showing Si-O pair-potential and its electrostatic and short-range components. A small separations it becomes attractive which is unphysical.

The first step to using ``exp_spline`` is to choose appropriate detachment and attachment points. This is perhaps best done plotting the two potential functions to be splined. The :ref:`as.buck <potable-buck>` and :ref:`as.zbl <potable-zbl` curves in :numref:`fig_spline_components`\ a, show that detachment and attachment at separations of 0.8 and 1.4 may be appropriate. 

.. figure:: figures/exp_spline_components.*
	:name: fig_spline_components
	:width: 60%
	:align: center

	The result of splining using the ``exp_spline`` potential.	

We now have enough information to spline the two potentials together. This gives is defined in the ``[Pair]`` section:

.. literalinclude:: example_files/exp_spline.aspot
	:lines: 6-12


The result of splining can be seen in :numref:`fig_spline_components`\ b (here it plotted with the Coulomb interaction included).

The complete ``potable`` input can be downloaded here: :download:`exp_spline.aspot <example_files/exp_spline.aspot>`\ .


.. _spline-buck4:

Buckingham-4 Spline ``buck4_spline``
------------------------------------

When used with ``buck4_spline`` the overall potential has the form:

.. math::
    :label: eq_buck4_spline

    V(r_{ij}) = 
    \begin{cases}
      V_\text{PotA}(r_{ij})                                                  & : 0 \leq r_{ij} \leq  r_\text{detach}\\
      a_0 + a_1 r_{ij} +a_2 r_{ij}^2+a_3 r_{ij}^3+a_4 r_{ij}^4+a_5 r_{ij}^5 & : r_\text{detach} < r_{ij} < r_\text{min}\\
      b_0 +b_1 r_{ij}+b_2 r_{ij}^2+b_3 r_{ij}^3                            & : r_\text{min} \leq r_{ij} < r_\text{attach}\\
      V_\text{PotB}(r_{ij})                                                   & : r_{ij} \geq r_\text{attach}\\
    \end{cases}


The spline form can be seen to be the combination of a fifth order and cubic polynomial. As for other spline types the  coefficients (:math:`a_{0 \ldots 5}` and :math:`b_{0 \ldots 3}`\ ) are chosen so that first and second derivatives (with respect to separation) are continuous through the attachment and detachment points. In addition, the coefficients are determined to give a stationary point where the two polynomials meet at :math:`r_\text{min}`\ .

The ``SPLINE_DEFN`` part of the ``spline()`` definition for ``buck4_spline`` is:

	| ``buck4_spline`` :math:`r_\text{min}`


Relationship to ``as.buck4`` potential-form
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As its name suggests, the ``buck4_spline`` is closely related to the :ref:`potform-buck4` potential-form. When used as a potential-form, the :math:`V_\text{PotA}(r_{ij})` and  :math:`V_\text{PotB}(r_{ij})` terms from eqn. :eq:`eq_buck4_spline` are pre-defined as the first and last terms of the :ref:`Buckingham <potable-buck>` potential:

.. math::

	V_\text{PotA}(r_{ij}) & = A_{ij} \exp\left(- \frac{r_{ij}}{\rho} \right) \\
	V_\text{PotB}(r_{ij}) & = - \frac{C}{r_{ij}} \\


The potential-form is therefore short-hand for the following::

	spline(as.buck A RHO 0.0 >R_DETACH buck4_spline R_MIN >R_ATTACH as.buck 0 1.0 C)

And would be specified as::

	as.buck4 A RHO C R_DETACH R_MIN R_ATTACH

Where: ``A``, ``RHO`` and ``C`` are the potential parameters and ``R_DETACH``, ``R_MIN`` and ``R_ATTACH`` define the splined region.

	
In most cases users will prefer the ``as.buck4`` form unless they have a specific reason to change the initial and final potential functions.


.. _spline-buck4_spline-example:

Example: redefining the Morelon model using ``buck4_spline``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this example the Morelon potential model will be re-implemented first using the ``spline()`` modifier then with the :ref:`as.buck4 <potform-buck4>`\ .

This model was discussed in a previous example: :ref:`potable-published-spline-coefficients`  where O-O and O-U interactions were defined like this:

.. literalinclude:: example_files/morelon.aspot
	:lines: 6-14


From this it can be seen that the O-U interaction is relatively simple whilst the O-O interaction four distinct parts (also defined in eqn. :eq:`eq_morelon-oo`\ ), and when examined more closely it can be seen it matches the four range Buckingham form.

The polynomials terms, with their pre-supplied coefficients, can be replaced by using the ``spline()`` modifier with the ``buck4_spline``\ spline type. This is now shown and can be downloaded as: :download:`morelon_buck4_spline.aspot <example_files/morelon_buck4_spline.aspot>`

.. literalinclude:: example_files/morelon_buck4_spline.aspot
	:lines: 6-13
	:linenos:


To change the starting potential and end potential lines 4 and 8 would be edited. To edit the detach and attach points then lines 5 and 7 would be altered. To change the position of the stationary point in the splined region, the value of 2.1 at the end of line 6 would be edited.

This input file could also be re-written to use the ``buck4`` potential-form, like this:

.. literalinclude:: example_files/morelon_buck4.aspot
	:lines: 6-8

This version of the model can be downloaded here: :download:`morelon_buck4.aspot <example_files/morelon_buck4.aspot>`

.. _pair_style soft:  https://lammps.sandia.gov/doc/pair_soft.html
.. _exprtk: http://www.partow.net/programming/exprtk/index.html