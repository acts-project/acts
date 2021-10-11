Python bindings for the Examples
================================

The examples part of ACTS ships with python bindings using the ``pybind11``
library. Building these bindings can be enabled via
``-DACTS_BUILD_EXAMPLES_PYTHON_BINDINGS=ON``, and requires a python installation
including the development files to be installed. You can then build the special
target `ActsPythonBindings` to build everything that can be accessed in python.

The build will create a setup script in ``$BUILD_DIR/python/setup.sh`` which
modifies ``$PYTHONPATH`` so that you can import the ``acts`` module in python.

Here is a minimal example of a python script using the example bindings, which
sets up the particle propagation and runs a few events.

.. code-block:: python

   import acts
   import acts.examples

   detector, trackingGeometry, contextDecorators  = acts.examples.GenericDetector.create()
   s = acts.examples.Sequencer(events=10)

   rnd = acts.examples.RandomNumbers(seed=42)

   nav = acts.Navigator(trackingGeometry=trackingGeometry)

   stepper = acts.EigenStepper(field)

   prop = acts.examples.ConcretePropagator(acts.Propagator(stepper, nav))

   alg = acts.examples.PropagationAlgorithm(
       propagatorImpl=prop,
       level=acts.logging.INFO,
       randomNumberSvc=rnd,
       ntests=1000,
       sterileLogger=True,
       propagationStepCollection="propagation-steps",
   )

   s.addAlgorithm(alg)

   s.addWriter(
       acts.examples.ObjPropagationStepsWriter(
           level=acts.logging.INFO,
           collection="propagation-steps",
           outputDir=outputDir + "/obj",
       )
   )

   s.addWriter(
       acts.examples.RootPropagationStepsWriter(
           level=acts.logging.INFO,
           collection="propagation-steps",
           filePath=outputDir + "/propagation_steps.root",
       )
   )

   s.run()

Python based example scripts
----------------------------

The repository contains a set of example scripts that can be used to execute various workflows.
They can be found in ``$REPO_ROOT/Examples/Scripts/Python``. Make sure you have run

.. code-block:: console

   source $BUILD_DIR/python/source.sh

to make sure python can find the ``acts`` module.

Python based unit tests
-----------------------

A number of unit tests based on the ``pytest`` library are shipped with the
repository. They are located under ``$REPO_ROOT/Examples/Python/tests``, and
intend to cover the public API of the python bindings. A set of tests also
executed the standalone example scripts.

To run these python based tests, ``pytest`` needs to be installed. It can be
installed via ``pip install pytest``. It is recommended to install this package
in `virtual environment`_. You can then simply run ``pytest`` from the
repository root.

.. _virtual environment: https://realpython.com/python-virtual-environments-a-primer/
