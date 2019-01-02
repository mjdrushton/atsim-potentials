.. _list-of-potential-forms:

***********************
List of Potential Forms
***********************

.. contents::
	:local:


.. rubric:: Key to features:

The **Features** field in the following potential description may contain the following values:

* **potential-form** - can be used in the ``[Pair]`` section of a :ref:`potable <potable-tool>` input file.
* **potential-function** - can also be used as a function in sections of input files that accept mathematical expressions.
* **deriv** - potential form provides an analytical derivative with respect to separation.
* **deriv2** - provides an analytical second derivative with respect to separation.


.. _potform-bornmayer:

Born-Mayer (bornmayer)
======================

.. math ::

	V(r_{ij}) = A \exp \left( \frac{-r_{ij}}{\rho} \right)


:potable signature: `as.bornmayer` :math:`A` :math:`\rho`
:Features: potential-form, potential-function, deriv, deriv2

.. seealso::

	* :ref:`potform-buck`
	


.. _potform-buck:

Buckingham (buck)
=================

Potential form due to R.A. Buckingham [Buckingham1938]_


.. math ::

	V(r_{ij}) = A \exp \left( - \frac{r_{ij}}{\rho} \right) - \frac{C}{r_{ij}^6}


:potable signature: ``as.buck`` :math:`A` :math:`\rho` :math:`C`
:Features: potential-form, potential-function, deriv, deriv2.

.. seealso::

	* :ref:`potform-bornmayer`
	* `Wikipedia - Buckingham potential <https://en.wikipedia.org/wiki/Buckingham_potential>`_


.. _potform-buck4:

Buckingham-4 (buck4)
====================

Four-range Buckingham potential due to B. Vessal et al. [Vessal1989]_\, [Vessal1993]_\ .

.. math::

    V(r_{ij}) = 
    \begin{cases}
      A \exp(-r_{ij}/\rho)                                                 , & 0 \leq r_{ij} \leq  r_\text{detach}\\
      a_0 + a_1 r_{ij} +a_2 r_{ij}^2+a_3 r_{ij}^3+a_4 r_{ij}^4+a_5 r_{ij}^5, & r_\text{detach} < r_{ij} < r_\text{min}\\
      b_0 +b_1 r_{ij}+b_2 r_{ij}^2+b_3 r_{ij}^3                           , & r_\text{min} \leq r_{ij} < r_\text{attach}\\
      -\frac{C}{r_{ij}^6}                                                  , & r_{ij} \geq r_\text{attach}\\
    \end{cases}

In other words this is a Buckingham potential in which the :ref:`Born-Mayer <potform-bornmayer>` component acts at small separations and the disprsion
term acts at larger separation. These two parts are linked by a fifth then third order polynomial (with a minimum formed in the spline
at :math:`r_\text{min}`).

The spline parameters (:math:`a_{0\ldots 5}` and :math:`b_{0\ldots 3}`) are subject to the constraints that :math:`V(r_{ij})`, first and second derivatives must be equal at the boundary points
and the function must have a stationary point at :math:`r_{min}`\ . The spline coefficients are automatically calculated by this potential-form.

.. note::

	Due to the complexity of calculating the spline-coefficients this potential form does not have an equivalent in the atsim.potentials.potentialfunctions
	module.

:potable signature: ``as.buck4`` :math:`A` :math:`\rho` :math:`C` :math:`r_\text{detach}` :math:`r_\text{min}` :math:`r_\text{attach}`
:Features: potential-form, deriv, deriv2

.. seealso::

		* :ref:`potform-buck4`
		* :ref:`aspot-splining`

.. _potform-constant:

Constant (constant)
===================

Potential form that always evaluates to a constant value.

.. math::

	V(r_{ij}) = C

:potable signature: ``as.constant`` C
:Features: potential-form, potential-function, deriv, deriv2

.. _potform-coul:

Coulomb (coul)
==============

Electrostatic interaction between two point charges.

.. math ::

      V(r_{ij}) = \frac{ q_i q_j }{4 \pi \epsilon_0 r_{ij} }

.. note:: 

	Constant value appropriate for :math:`r_{ij}` in angstroms and energy in eV.


:potable signature: ``as.coul`` :math:`q_i` :math:`q_j`
:Features: potential-form, potential-function, deriv, deriv2

.. _potform-exp_spline:

Exponential Spline (exp_spline)
===============================

Exponential spline function (as used in splining routines).

      .. math::

              V(r_{ij}) = \exp \left( B_0 + B_1 r_{ij} + B_2 r_{ij}^2 + B_3 r_{ij}^3 + B_4 r_{ij}^4 + B_5 r_{ij}^5 \right) + C

Where :math:`B_0`\ , :math:`B_1`\ , :math:`B_2`\ , :math:`B_3`\ , :math:`B_4`\ , :math:`B_5`\ , :math:`C` are spline coefficients.

:potable signature: ``as.exp_spline`` :math:`B_0` :math:`B_1` :math:`B_2` :math:`B_3` :math:`B_4` :math:`B_5` :math:`C`
:Features: potential-form, potential-function, deriv, deriv2

.. seealso::

	* :ref:`aspot-splining`



.. _potform-hbnd:

Hydrogen Bond 12-10 (hbnd)
==========================

.. math::

      V(r_{ij}) = \frac{A}{r_{ij}^{12}} - \frac{B}{r_{ij}^{10}}

:portable signature: ``as.hbnd`` :math:`A` :math:`B`
:Features: potential-form, potential-function, deriv, deriv2

.. _potform-lj:

Lennard-Jones 12-6 (lj)
=======================

Potential form first proposed by John Lennard-Jones in 1924 [Lennard-Jones1924]_\ .

.. math ::

	V(r_{ij}) = 4 \epsilon \left( \frac{\sigma^{12}}{r_{ij}^{12}} - \frac{\sigma^6}{r_{ij}^6} \right)

:math:`\epsilon` defines depth of potential well and :math:`\sigma` is the separation at which :math:`V(r_{ij})` is zero.

:potable signature: ``as.lj`` :math:`\epsilon` :math:`\sigma`
:Features: potential-form, potential-function, deriv, deriv2

.. seealso::

	* `Wikipedia - Lennard-Jones potential <https://en.wikipedia.org/wiki/Lennard-Jones_potential>`_


.. _potform-morse:

Morse (morse)
=============

.. math ::

	V(r_{ij}) = D  \left[ \exp \left( -2 \gamma (r_{ij} - r_*) \right) - 2 \exp \left( -\gamma (r - r_*) \right) \right]

:math:`-D` is the potential well depth at an equilibrium separation of :math:`r_*`\ .

:potable signature: ``as.morse`` :math:`\gamma` :math:`r_*` :math:`D`
:Features: potential-form, potential-function, deriv, deriv2

.. seealso::

	* `Wikipedia - Morse potential <https://en.wikipedia.org/wiki/Morse_potential>`_



.. _potform-polynomial:

Polynomial (polynomial)
=======================

Polynomial of arbitrary order.

.. math::

	V(r_{ij}) = C_0 + C_1 r_{ij} + C_2 r_{ij}^2 + \dots + C_n r_{ij}^n

This function accepts a variable number of arguments which are :math:`C_0, C_1, \dots, C_n` respectively.

:potable signatures: ``as.polynomial`` :math:`C_0 ... C_n`
:Features: potential-form, potential-function, deriv, deriv2


.. _potform-sqrt:

Square Root (sqrt)
==================

Potential function is:

.. math::

	U(r_{ij}) = G\sqrt{r_{ij}}

:potable signature: ``as.sqrt`` :math:`G`
:Features: potential-form, potential-function, deriv, deriv2


.. _potform-tang-toennies:

Tang-Toennies (tang_toennies)
=============================

This potential form was derived to describe the Van der Waal's interactions between the noble gases (He to Rn) by Tang and Toennies [Tang2003]_\ .

This has the following form:

.. math::

	V(r) = A \exp(-br) - \sum_{n=3}^N f_{2N} (bR) \frac{C_{2N}}{R^{2N}}

Where:

.. math::

	f_{2N}(x) = 1- \exp(-x) \sum_{k=0}^{2n} \frac{x^k}{k!}

:potable signature: ``as.tang_toennies`` :math:`A` :math:`b` :math:`C_6` :math:`C_8` :math:`C_{10}`
:Features: potential-form, potential-function, deriv, deriv2


.. _potform-zero:

Zero (zero)
===========

Potential form which returns zero for all separations.

.. math::

	V(r) = 0

:potable signature: ``as.zero``
:Features: potential-form, potential-function, deriv, deriv2


.. _potform-zbl:

Ziegler-Biersack-Littmark (zbl)
===============================

Ziegler-Biersack-Littmark screened nuclear repulsion for describing high energy interactions [Ziegler2015]_\ .

.. math::

	V(r)    & = \frac{1}{4\pi\epsilon_0} \frac{Z_1}{Z_2} \phi(r/a) + S(r) \\
	a       & = \frac{0.46850}{Z_i^{0.23} + Z_j^{0.23}} \\
	\phi(x) & = 0.18175 \exp(-3.19980x) + 0.50986 \exp(-0.94229x) + 0.28022\exp(-0.40290x) + 0.02817\exp(-0.20162x) \\


Where :math:`Z_i` and :math:`Z_j` are the atomic numbers of two species.

:potable signature: ``as.zbl`` :math:`Z_i` :math:`Z_j`
:Features: potential-form, potential-function, deriv, deriv2

