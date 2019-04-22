
.. _ref-potable-input-format:

************************
``potable`` input format
************************

.. _ref-potable-input-tabulation:

[Tabulation]
============

.. _ref-potable-input-tabulation-cutoff:

cutoff
------

:Item: ``cutoff``
:Format: float
:Description: Defines upper bound of functions tabulated in terms of separation. This is used in a pair with ns tabulated in terms of separation. This directive is used together with :ref:`nr (number of rows) <ref-potable-input-tabulation-nr>` or :ref:`dr  (step size) <ref-potable-input-tabulation-dr>` to give the extent and resolution of a tabulated function.


.. _ref-potable-input-tabulation-dr:

dr
--

:Item: ``dr``
:Format: float
:Description: Defines the step size between rows of functions tabulated in terms of separation. This directive is used together with :ref:`nr <ref-potable-input-tabulation-nr>` or :ref:`cutoff <ref-potable-input-tabulation-cutoff>` to define resolution and extent of these functions.


.. _ref-potable-input-tabulation-nr:

nr
--

:Item: ``nr``
:Format: int
:Description: Defines the number of rows when functions are tabulated in terms of separation. This directive is used either with :ref:`dr <ref-potable-input-tabulation-dr>` or :ref:`cutoff <ref-potable-input-tabulation-cutoff>` to give the range and resolution of the tabulated function.

.. _ref-potable-input-tabulation-target:

target
------

:Item: ``target``
:Format: str
:Valid Options: ``DL_POLY|DLPOLY``, ``DLPOLY_EAM_fs`` ``DLPOLY_EAM``, ``GULP``, ``LAMMPS_eam_alloy|setfl``, ``LAMMPS``, ``setfl_fs``
:Description: Specifies the format that tabulation will be written in.


.. _ref-potable-input-pair:

[Pair]
======

Pair-potentials are defined in this section of the file. See :ref:`potable-pair-section` for full description.

.. seealso::
  
  See also:

  * Potential-forms are parametrised in this section:
      - :ref:`list-of-potential-forms` - reference list of pre-defined potential forms.
      - Custom functions are defined in the ``[Potential-Form]`` section:

          + :ref:`potable-potential-form` - custom potential-forms are introduced here.
          + :ref:`ref-potential-form` - reference information for ``[Potential-Form]`` section.
 
  * Potential-modifiers may be used in this section:

      - :ref:`potential-modifiers` - are introduced here.
      - :ref:`ref-potential-modifiers` - list of potential-modifiers.



[Potential-Form]
================

Custom functional forms are defined in this section. See :ref:`potable-potential-form` where it is introduced.

.. seealso::

    * The syntax used by the mathematical expressions defined in the ``[Potential-Form]`` is `defined here <http://www.partow.net/programming/exprtk/index.html>`_\ .

.. _ref-potable-input-pymath:

Python maths functions supported in mathematical expressions
------------------------------------------------------------

The mathematical expressions used in the ``[Potential-Form]`` section of ``potable`` input allow a subset of functions from the `math <https://docs.python.org/3/library/math.html>`_ module to be used. These are accesible via the ``pymath.*`` namespace prefix. An example of this is provided here: :ref:`potable-potential-form-formula-syntax`

The list of functions accessible through ``pymath.*`` are below. In general, functions that return multiple values do not appear:

    * `acos(x) <https://docs.python.org/3/library/math.html#math.acos>`_
    * `acosh(x) <https://docs.python.org/3/library/math.html#math.acosh>`_
    * `asinh(x) <https://docs.python.org/3/library/math.html#math.asinh>`_
    * `atan(x) <https://docs.python.org/3/library/math.html#math.atan>`_
    * `atan2(x,y) <https://docs.python.org/3/library/math.html#math.atan2>`_
    * `atanh(x) <https://docs.python.org/3/library/math.html#math.atanh>`_
    * `cos(x) <https://docs.python.org/3/library/math.html#math.cos>`_
    * `cosh(x) <https://docs.python.org/3/library/math.html#math.cosh>`_
    * `degrees(x) <https://docs.python.org/3/library/math.html#math.degrees>`_
    * `exp(x) <https://docs.python.org/3/library/math.html#math.exp>`_
    * `factorial(x) <https://docs.python.org/3/library/math.html#math.factorial>`_
    * `fsum(*args) <https://docs.python.org/3/library/math.html#math.fsum>`_ 

        + This function is called slightly differently than in native Python.
        + In Python you pass in a single iterable to this function. This expression: ``math.fsum([1,2,3,4])`` would be written ``pymath.fsum(1,2,3,4)`` in a ``potable`` formula.

    * `gcd(a,b) <https://docs.python.org/3/library/math.html#math.gcd>`_
    * `hypot(x,y) <https://docs.python.org/3/library/math.html#math.hypot>`_
    * `ldexp(a,b) <https://docs.python.org/3/library/math.html#math.ldexp>`_
    * `log(*args) <https://docs.python.org/3/library/math.html#math.log>`_
    * `log10(x) <https://docs.python.org/3/library/math.html#math.log10>`_
    * `log1p(x) <https://docs.python.org/3/library/math.html#math.log1p>`_
    * `log2(x) <https://docs.python.org/3/library/math.html#math.log2>`_
    * `pow(x,a) <https://docs.python.org/3/library/math.html#math.pow>`_
    * `radians(x) <https://docs.python.org/3/library/math.html#math.radians>`_
    * `sin(x) <https://docs.python.org/3/library/math.html#math.sin>`_
    * `sinh(x) <https://docs.python.org/3/library/math.html#math.sinh>`_
    * `sqrt(x) <https://docs.python.org/3/library/math.html#math.sqrt>`_
    * `sqrt(x) <https://docs.python.org/3/library/math.html#math.sqrt>`_
    * `tan(x) <https://docs.python.org/3/library/math.html#math.tan>`_
    * `tanh(x) <https://docs.python.org/3/library/math.html#math.tanh>`_
    * `trunc(x) <https://docs.python.org/3/library/math.html#math.trunc>`_