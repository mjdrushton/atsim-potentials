
****************************************
Working with ``potable`` files in Python
****************************************

Using :class:`~atsim.potentials.config.Configuration` class to create ``Pair_Tabulation`` and ``EAM_Tabulation`` objects from ``potable`` input
===============================================================================================================================================

:class:`atsim.potentials.config.Configuration` is a factory class which accepts ``potable`` input and uses it to create tabulation objects.

Tabulation objects are typically created by passing a file like object containing ``potable`` input to the :meth:`~atsim.potentials.config.Configuration.read` method of :class:`~atsim.potentials.config.Configuration`\ . 


.. _python_potable_api_example:

Example
-------

The following example demonstrates how to create a :class:`atsim.potentials.pair_tabulation.Pair_Tabulation` object from ``potable`` input by using the :class:`atsim.potentials.config.Configuration` class.

The aim of the example is to find the spline coefficients for following potable input (described elsewhere :ref:`aspot-splining`\ )

.. literalinclude:: example_files/exp_spline.aspot



Running the following script :download:`python_potable_api.py <example_files/python_potable_api.py>` will print out the spline coefficients:

.. literalinclude:: example_files/python_potable_api.py




Overriding and adding items
===========================

Items can be amended or added to the ``potable`` input before it is passed to the :class:`~atsim.potentials.config.Configuration` class. This is done by passing a :class:`atsim.potentials.config.ConfigParser` to the :meth:`~atsim.potentials.config.Configuration.read_from_parser` method of the :class:`~atsim.potentials.config.Configuration` object. 

The constructor of :class:`~atsim.potentials.config.ConfigParser` accepts ``overrides`` and ``additional`` parameters, each of which accept lists of :class:`atsim.potentials.config.ConfigParserOverrideTuple`\ .

:class:`~atsim.potentials.config.ConfigParserOverrideTuple` is a :class:`collections.namedtuple` with three properties ``section``\ , ``key`` and ``value``\ . The first two uniquely identify a location in the ``potable`` input whilst ``value`` specifies what should be added or changed.

So to add an additional pair-potential to potable input contained in a file given by ``fp`` the :class:`~atsim.potentials.config.ConfigParser` would be defined as:

.. code::

    cp = ConfigParser(fp,
                      additional=[
                          ConfigParserOverrideTuple(
                              "Pair", "O-O", "as.buck 444.7686 0.402 0.0")
                      ])


This would be the same as if the following had been given in the original ``potable`` input::

    [Pair]
    O-O : as.buck 444.7686 0.402 0.0



Similarly to change the tabulation target of a potable file you could use:

.. code::

    cp = ConfigParser(fp,
                      overrides=[
                          ConfigParserOverrideTuple(
                              "Tabulation", "target", "LAMMPS"),
                      ])


A tabulation object is then obtained by combining the :class:`~atsim.potentials.config.ConfigParser` with :class:`~atsim.potentials.config.Configuration`\ :

.. code::

    tabulation = Configuration().read_from_parser(cp)


.. _python_potable_api_override_example:

Example
-------

This example shows the use of overrides and additional items. Again the potable input from :ref:`aspot-splining` is used.

.. literalinclude:: example_files/exp_spline.aspot


In :download:`python_potable_api.py <example_files/python_potable_api.py>` the tabulation target and potential cutoff in the ``[Tabulation]`` section are overriden. An additional ``O-O`` interaction is added to the ``[Pair]`` section. 

This is then used to create a tabulation object which is finally output to the screen:

.. literalinclude:: example_files/python_potable_override.py





