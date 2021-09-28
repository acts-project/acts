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