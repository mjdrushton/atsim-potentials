
.. _ref-potable-input-format:

************************
``potable`` input format
************************


.. _ref-potable-eam-density:

[EAM-Density]
=============

The density functions for embedded atom models are specified in this section. The input takes different forms depending on whether the standard embedded atom model or Finnis-Sinclair variant are being used. 

Both standards have the following general form::

    INTERACTION : POTENTIAL_FORM PARAM_1 PARAM_2 ... PARAM_N

Where:

    * ``POTENTIAL_FORM PARAM_1 ...`` \ : density functions use the same rules to instantiate potential forms as in the :ref:`[Pair] <ref-potable-input-pair>` section.

    * ``INTERACTION`` specifies the density this potential-form represents:

        * **Standard EAM:** standard EAM uses the same  function for the density surrounding any central atom of any given species. Consequently in these cases ``INTERACTION`` is a single species label. So the density function of aluminium would take the form::
        
            [EAM-Density]
            Al : POTENIAL_FORM ...


        * **Finnis-Sinclair:** in this variant of EAM density functions change depending on types of the cental atom and surrounding atom to be embedded. Consequently the following form is used::
        
            [EAM-Density]
            A->B : POTENTIAL_FORM ...


         * Where ``A`` is the central atom type and ``B`` is the type of the embedding atom. To define a density function for the density of nickel being embedded at an aluminium site this would be used::
        
            [EAM-Density]
            Al->Ni : POTENTIAL_FORM ...

        
         * It should be noted that ``A->B`` and ``B->A`` must be specified separately even if the same density function is used for both. If not given null (i.e. ``as.zero``\ ) density functions are implicitly defined for missing interactions.

.. seealso::

    * _potable-many-body-models


.. _ref-potable-eam-embed:

[EAM-Embed]
===========

Embedding functions for many-body models are defined in this section. 

Entries have the following form::

    SPECIES : POTENTIAL_FORM PARAM_1 PARAM_2 ... PARAM_N

Where:

    * ``SPECIES`` is atomic type at which the surrounding electron density will be embedded using the specified potential form.
    * ``POTENTIAL_FORM PARAM_1 ...`` \ : embedding functions instantiate potential forms in the same way as in the :ref:`[Pair] <ref-potable-input-pair>` section.


.. note::

    Embedding functions are tabulated using rho values. The resolution and extent of functions in rho are defined by ``drho``\ , ``nrho`` and ``cutoff_rho`` in the :ref:`[Tabulation] section <ref-potable-input-tabulation>`\ .


.. seealso::

    * _potable-many-body-models


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
          + :ref:`ref-potable-potential-form` - reference information for ``[Potential-Form]`` section.
 
  * **Potential-modifiers are described in thiese sections:**

      - :ref:`potential-modifiers` - are introduced here.
      - :ref:`list-of-potential-modifiers` - list of potential-modifiers.


.. _ref-potable-potential-form:

[Potential-Form]
================

Custom functional forms are defined in this section. See :ref:`potable-potential-form` where it is introduced.

.. seealso::

    * The syntax used by the mathematical expressions defined in the ``[Potential-Form]`` is `defined here <http://www.partow.net/programming/exprtk/index.html>`_\ .

.. _ref-potable-input-pymath:

Python maths functions supported in mathematical expressions
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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

.. _ref-potable-input-tabulation:

[Tabulation]
============

The section of the input file which defines how a model should be tabulated.

Fields
++++++

.. _ref-potable-input-tabulation-cutoff:

cutoff
------

:Item: ``cutoff``
:Format: float
:Description: Defines upper bound of functions tabulated in terms of separation. This is used in a pair with ns tabulated in terms of separation. This directive is used together with :ref:`nr (number of rows) <ref-potable-input-tabulation-nr>` or :ref:`dr  (step size) <ref-potable-input-tabulation-dr>` to give the extent and resolution of a tabulated function.

.. _ref-potable-input-tabulation-cutoff-rho:

cutoff_rho
----------

:Item: ``cutoff_rho``
:Format: float
:Description: Used to define cutoff for functions tabulated in terms of electron density (rho) e.g. for :ref:`ref-potable-eam-embed` functions. This option defines the upper bound of rho values included in the tabulation of these functions. This directive is used together with :ref:`nrho <ref-potable-input-tabulation-nrho>` or :ref:`cutoff_rho <ref-potable-input-tabulation-cutoff-rho>` to define resolution and extent of density functions.


.. _ref-potable-input-tabulation-dr:

dr
--

:Item: ``dr``
:Format: float
:Description: Defines the step size between rows of functions tabulated in terms of separation. This directive is used together with :ref:`nr <ref-potable-input-tabulation-nr>` or :ref:`cutoff <ref-potable-input-tabulation-cutoff>` to define resolution and extent of these functions.


.. _ref-potable-input-tabulation-drho:

drho
----

:Item: ``drho``
:Format: float
:Description: Used to define resolution of functions tabulated in terms of electron density (rho) e.g. for :ref:`ref-potable-eam-embed` functions. This option defines the rho increment for such functions. This directive is used together with :ref:`nrho <ref-potable-input-tabulation-nrho>` or :ref:`cutoff_rho <ref-potable-input-tabulation-cutoff-rho>` to define resolution and extent of these functions.


.. _ref-potable-input-tabulation-nr:

nr
--

:Item: ``nr``
:Format: int
:Description: Defines the number of rows when functions are tabulated in terms of separation. This directive is used either with :ref:`dr <ref-potable-input-tabulation-dr>` or :ref:`cutoff <ref-potable-input-tabulation-cutoff>` to give the range and resolution of the tabulated function.


.. _ref-potable-input-tabulation-nrho:

nrho
----

:Item: ``nrho``
:Format: int
:Description: Used to define cutoff (in conjunction with ``drho``\ ) for functions tabulated in terms of electron density (rho) e.g. for :ref:`ref-potable-eam-embed` functions. This option defines the number of rho values included in the tabulation of these functions. This directive is used together with :ref:`nrho <ref-potable-input-tabulation-nrho>` or :ref:`cutoff_rho <ref-potable-input-tabulation-cutoff-rho>` to define resolution and extent of density functions.


.. _ref-potable-input-tabulation-target:

target
------

:Item: ``target``
:Format: str
:Valid Options: ``DL_POLY|DLPOLY``, 
    ``DLPOLY_EAM_fs``,
    ``DLPOLY_EAM``, 
    ``excel``,
    ``excel_eam``,
    ``excel_eam_fs``,
    ``GULP``, 
    ``LAMMPS_eam_alloy|setfl``,
    ``LAMMPS``, 
    ``setfl_fs``
:Description: Specifies the format that tabulation will be written in.


[Table-Form]
============

The ``[Table-Form]`` section is used to define functions from pre-tabulated data that may be used in the same way as a custom ``[Potentia-Form]``\ . Data is specified using the ``x`` and ``y`` options or the ``xy`` option.

To provide a continuous function interpolation is performed between data points, the interpolation method is set using the ``interpolation`` option.

Naming Table Form
+++++++++++++++++

To allow a ``[Table-Form]`` to be used in sections such as ``[Pair]``\ , ``[EAM-Embed]`` and ``[EAM-Density]`` it is necessary to give it a unique label. This is done by including it in the section header following a colon::

    [Table-Form:NAME]

Therefore to create a ``[Table-Form]`` named ``tabulated`` the following definition could be used::

    [Table-Form:tabulated]
    interpolation: cubic_spline
    x : 0.0 1.0 2.0 3.0
    y : 0.0 2.0 3.0 4.0


This could then be referenced in another section using this name. e.g. ::

    [Pair]
    Si-O : tabulated


Fields
++++++

.. _ref-potable-input-table-form-interpolation:

interpolation
-------------

:Item: ``interpolation``
:Format: Currently this option only accepts ``cubic_spline``
:Description: Sets interpolation type.


.. _ref-potable-input-table-form-x:

x
-

:Item: ``x``
:Format: List of space separated float values.
:Description: Define x values of tabulated data. Must be used with ``y`` option.
:Example: To define a linear function the following could be used:

::

    [Table-Form:linear]
    interpolation: cubic_spline
    x : 0.0 1.0 2.0 3.0
    y : 0.0 2.0 3.0 4.0


.. _ref-potable-input-table-form-xy:

xy
--

:Item: ``xy``
:Format: List of space separated float values.
:Description: Allows x and y values of data to be specified as series of pairs.
:Example: To define a linear function the following could be used:

::

    [Table-Form:linear]
    interpolation: cubic_spline
    xy: 0.0 0.0
        1.0 2.0
        2.0 3.0
        3.0 4.0

y
-

:Item: ``y``
:Format: List of space separated float values.
:Description: Define y values of tabulated data. Must be used with ``x`` option.
:Example: See documentation for :ref:`ref-potable-input-table-form-x` option.
