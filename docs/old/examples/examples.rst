Examples
========

ACTS ships with a comprehensive set of examples. These examples leverage a
custom event-processing framework, that is expressly **not intended to be used
in any kind of production environment**. These examples demonstrate how to set
up and configure different components in the ACTS library to assemble a track
reconstruction chain.


At the time of writing, there are two aspects to the ACTS examples:

#. Example executables for different purposes. This is the original form of
   examples provided by ACTS. A large number of separate executables are be built
   if ``-DACTS_BUILD_EXAMPLES=ON``, the exact set is also influenced by which
   plugins are enabled in the build. These executables are configured by a number of
   command line options, for example to set the number of events to be processed,
   or which output formats to read from / write to.

#. Standalone Performance and Analysis applications based on ROOT. These applications
   are built on top of the ROOT based output writers of the ``Examples`` folder, they
   comprise of track reconstruction performance validation and material validation
   code.

#. Python bindings for the various components of the examples. These bindings
   allow more flexible combination of the example components into scripts, and can
   overcome the complexity of the large number of executables and command line
   options. The idea is that these scripts will serve as true examples, where
   modifications to the actual python code will be easy, and encouraged.


.. note:: This section of the documentation contains a set of *how-to* guides,
          which describe different example executables, and how to combine them
          with one another to assemble several workflows.

.. warning:: Please, do not use the event-processing framework or the example
             executables in a production environment! They are meant for demonstration
             and prototyping purposes only!

.. toctree::
   :maxdepth: 2

   Python bindings <python_bindings>
   OpenDataDetector full chain <full_chain_odd>
   alignment
   How-to guides <howto/howto>
