**********
Python API
**********

.. todo::

	Update this to reflect new API structure


.. todo::

	Potential object doesn't talk about handling of derivatives correctly - this needs to be updated.




.. py:module:: atsim.potentials

``atsim.potentials``
====================

Pair Potential Tabulation
-------------------------

.. autoclass:: atsim.potentials.Potential
    :members:
    :undoc-members:

    .. automethod:: __init__
    

.. autofunction:: atsim.potentials.writePotentials


Embedded Atom Method Tabulation
-------------------------------

.. autoclass:: atsim.potentials.EAMPotential
    :members:
    :undoc-members:

    .. automethod:: __init__

.. autofunction:: atsim.potentials.writeFuncFL

.. autofunction:: atsim.potentials.writeSetFL

.. autofunction:: atsim.potentials.writeSetFLFinnisSinclair

.. autofunction:: atsim.potentials.writeTABEAM 

.. autofunction:: atsim.potentials.writeTABEAMFinnisSinclair 


Miscellaneous Functions
-----------------------

.. autofunction:: atsim.potentials.plot

.. autofunction:: atsim.potentials.plotToFile

.. autofunction:: atsim.potentials.plotPotentialObject

.. autofunction:: atsim.potentials.plotPotentialObjectToFile

..  autofunction:: atsim.potentials.plus

.. autoclass:: atsim.potentials.SplinePotential
    :members:
    :undoc-members:

    .. automethod:: __init__


.. autoclass:: atsim.potentials.TableReader
    :members:
    :undoc-members:



.. _atsim_potentials_potentialforms:

``atsim.potentials.potentialforms``
===================================

.. automodule:: atsim.potentials.potentialforms
    :members:
    :undoc-members:

