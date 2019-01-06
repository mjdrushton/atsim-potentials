.. _list-of-potential-modifiers:

***************************
List of Potential Modifiers
***************************

Potential modifiers are described here: :ref:`potential-modifiers`\ .

A list of available potential modifiers is provided here.


.. _modifier-spline:

``spline()``
============

Modifier that smoothly splines between two potential forms by linking them with an intermediate spline.

``spline()`` takes a single argument which is  defined as a :ref:`multi-range potential <multi-range-potentials>`. This must define three ranges:

#. Start potential
#. Interpolating spline
#. End potential

The **Interpolating spline** section has the form::

	SPLINE_LABEL SPLINE_PARAMETERS

Where the ``SPLINE_LABEL`` defines the type of spline to be used and the (optional) ``SPLINE_PARAMETERS`` is a list of space separated options taken by the spline function.

A list of spline types usable with ``SPLINE_LABEL`` is now given:

``buck4_spline``
----------------

:Spline Signature: ``buck4_spline`` :math:`r_\text{min}`
:Description: Combination of a fifth and third order polynomial joined by a stationary point at :math:`r_\text{min}`.
	This is the spline used in the well-known :ref:`four-range Buckingham potential form <potform-buck4>`\ .
:See also: 
	* :ref:`spline-buck4`
	* :ref:`potform-buck4`

\

``exp_spline``
--------------

:Spline Signature: ``exp_spline``
:Description: Exponential of fifth order polynomial.
:See also: 
	* :ref:`spline-exp_spline`
	* :ref:`potform-exp_spline`


Example:
--------
  
A configuration string might be defined as::

	[Pair]

	Si-O : spline(>0 as.zbl 14 8 >=0.8 exp_spline >=1.4 as.buck 180003 0.3 32.0)


This would create a :ref:`zbl <potform-zbl>` and :ref:`Buckingham <potform-buck>` potential connected by an exponential spline when `r` is between 0.8 and 1.4.

.. seealso::

	* Splining is introduced in more detail here: :ref:`aspot-splining`\ .
	* List of examples:


	    - :ref:`spline-exp_spline-example`\ .
	    - :ref:`spline-buck4_spline-example`\ .


.. _modifier-sum:

``sum()``
=========

Modifier that sums all the potentials given as arguments. 

If the potentials provide analytical derivatives the ``sum()`` modifier will combine these correctly.

Example:
--------

Any number of potential instances can be summed::

	[Pair]
	# Evaluates to 3
	A-A : sum(as.constant 1.0, as.constant 2.0)
	# Evaluates to 6
	B-B : sum(as.constant 1.0, as.constant 2.0, as.constant 3.0)

.. seealso::

	* This modifier is used in the following examples:

	    - :ref:`quick-start`
	    - :ref:`potable-potential-form-basak-example`



