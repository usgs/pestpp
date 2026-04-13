Run Managers
============

The run manager layer handles model execution — serial, distributed (Panther/YAMR),
external, and language-wrapper modes.

Abstract Base
-------------

.. doxygenfile:: RunManagerAbstract.cpp
   :project: pestpp

.. doxygenfile:: RunManagerAbstract.h
   :project: pestpp

.. doxygenfile:: RunStorage.cpp
   :project: pestpp

.. doxygenfile:: RunStorage.h
   :project: pestpp

.. doxygenfile:: Serialization.cpp
   :project: pestpp

.. doxygenfile:: Serialization.h
   :project: pestpp

.. doxygenfile:: model_interface.cpp
   :project: pestpp

.. doxygenfile:: model_interface.h
   :project: pestpp

.. doxygenfile:: debug.cpp
   :project: pestpp

.. doxygenfile:: debug.h
   :project: pestpp

Serial
------

.. doxygenfile:: RunManagerSerial.cpp
   :project: pestpp

.. doxygenfile:: RunManagerSerial.h
   :project: pestpp

Panther (YAMR Distributed)
---------------------------

.. doxygenfile:: RunManagerPanther.cpp
   :project: pestpp

.. doxygenfile:: RunManagerPanther.h
   :project: pestpp

.. doxygenfile:: PantherAgent.cpp
   :project: pestpp

.. doxygenfile:: PantherAgent.h
   :project: pestpp

External
--------

.. doxygenfile:: RunManagerExternal.cpp
   :project: pestpp

.. doxygenfile:: RunManagerExternal.h
   :project: pestpp

Language Wrappers
-----------------

.. doxygenfile:: RunManagerCWrapper.cpp
   :project: pestpp

.. doxygenfile:: RunManagerCWrapper.h
   :project: pestpp

.. doxygenfile:: RunManagerFortranWrapper.cpp
   :project: pestpp

.. doxygenfile:: RunManagerFortranWrapper.h
   :project: pestpp
