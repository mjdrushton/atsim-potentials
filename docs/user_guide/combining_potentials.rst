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


Example: parametrising a model using published spline coefficients
------------------------------------------------------------------

The UO\ :sub:`2` potential model by Morelon describes U-O and O-O interactions [Morelon2003]_\ . The former employs the :ref:`potform-bornmayer` potential-form while the latter uses the combination of a 3rd and 5th order polynomial spline to link an :math:`A \exp\left( - \frac{r_{ij}}{\rho}\right)` term to a :math:`-\frac{C}{r_{ij}^6}` term. The spline-coefficients for the Morelon model are available in a paper by Potashnikov et al. [Potashnikov2011]_\ . They give the O-O interaction as:

.. math::

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

.. todo:

  Write section on splines


.. _pair_style soft:  https://lammps.sandia.gov/doc/pair_soft.html
.. _exprtk: http://www.partow.net/programming/exprtk/index.html