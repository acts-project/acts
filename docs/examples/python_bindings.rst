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

   import os

   import acts
   import acts.examples

   detector, trackingGeometry, contextDecorators  = acts.examples.GenericDetector.create()
   s = acts.examples.Sequencer(events=10)

   rnd = acts.examples.RandomNumbers(seed=42)

   nav = acts.Navigator(trackingGeometry=trackingGeometry)

   field = acts.ConstantBField(acts.Vector3(0, 0, 2 * acts.UnitConstants.T))
   stepper = acts.EigenStepper(field)

   prop = acts.examples.ConcretePropagator(acts.Propagator(stepper, nav))

   alg = acts.examples.PropagationAlgorithm(
       propagatorImpl=prop,
       level=acts.logging.INFO,
       randomNumberSvc=rnd,
       ntests=1000,
       sterileLogger=False,
       propagationStepCollection="propagation-steps",
   )

   s.addAlgorithm(alg)

   outputDir = "."
   objDir = outputDir + "/obj"
   if not os.path.exists(objDir):
      os.mkdir(objDir)
   
   s.addWriter(
       acts.examples.ObjPropagationStepsWriter(
           level=acts.logging.INFO,
           collection="propagation-steps",
           outputDir=objDir,
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

   source $BUILD_DIR/python/setup.sh

to make sure python can find the ``acts`` module.

Python based unit tests
-----------------------

A number of unit tests based on the ``pytest`` library are shipped with the
repository. They are located under ``$REPO_ROOT/Examples/Python/tests``, and
intend to cover the public API of the python bindings. A set of tests also
executed the standalone example scripts.

To run these python based tests, ``pytest`` and a few other dependencies need to be installed. They can be
installed via ``pip install -r Examples/Python/tests/requirements.txt`` from the repository root. It is recommended to install these packages
in `virtual environment`_. You can then simply run ``pytest`` from the
repository root.

.. _virtual environment: https://realpython.com/python-virtual-environments-a-primer/
