.. _list-of-potential-modifiers:

***************************
List of Potential Modifiers
***************************

Potential modifiers are described here: :ref:`potential-modifiers`\ .

A list of available potential modifiers is provided here.

.. _modifier-pow:

``pow()``
=========

Modifier that raises each potential-form to the power of the next.

If the potentials provide analytical derivatives, the ``pow()`` modifier will combine these correctly.

Example:
--------

To take the square of the :ref:`modifier-sum` of a series of potential forms you could use::

	[Pair]
	# This would evaluate to 0.16
	A-B : pow(sum(as.constant -1, as.constant 0.1, as.constant 0.5),
	          as.constant 2)


``pow()`` can take more than two potential forms as its arguments::

	[Pair]
	# This would evaluate to 2^(2^3) = 256 
	A-B : pow(as.constant 2, as.constant 3, as.constant 2)


You aren't restricted to using constant values as arguments::

	[Pair]
	# This is equivalent to 2^(0.5r + r^2)
	A-B : pow(as.constant 2, as.polynomial 0 0.5 1)

.. _modifier-product:

``product()``
=============

Modifier that takes the product of the potential-forms provided to it as arguments.

If the potentials provide analytical derivatives the ``product()`` modifier will combine these correctly.

Example:
--------

Any number of potential instances can be multiplied by each other::

	[Pair]
	# Evaluates to 16
	A-A : product(as.constant 2.0, as.constant 2.0, as.constant 4.0)
	
	# Apply a soft-cutoff at 2.5 Angs to a Buckingham potential
	# This defines a custom function in the [Potential-Form] section
	# based on the complementary error function for this purpose. 
	B-B : product(as.buck 1000.0 0.2 32.0,
	              truncate 2.5)

	[Potential-Form]
	truncate(rij, cutoff) = erfc(4*(rij-cutoff))/2.0

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



