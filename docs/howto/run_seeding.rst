Run the seeding example
===============================

Prerequisites
-------------

ACTS needs to be built with ``ACTS_BUILD_EXAMPLES=ON`` to enable this example.
Before running the seeding example, you have to prepare a simulated dataset. 
The event generation and fatras simulation are explained in :ref:`simulate-TrackML`.


Run the seeding example
----------------------


Suppose the generated ttbar sample is available at ``data/sim_generic/ttbar_pu200``.
The seed finding example reads the fatras simulation output and creates smeared measurements from the truth hits.
After that, it creates space points from the smeared measurements, and the seed finding algorithm is performed.
The following parameters for the seed finding are set in ``/acts/Examples/Run/Seeding/Common/SeedingExample.cpp``.

.. code-block:: console

  seedingCfg.rMax = 200.;
  seedingCfg.deltaRMax = 60.;
  seedingCfg.collisionRegionMin = -250;
  seedingCfg.collisionRegionMax = 250.;
  seedingCfg.zMin = -2000.;
  seedingCfg.zMax = 2000.;
  seedingCfg.maxSeedsPerSpM = 1;
  seedingCfg.cotThetaMax = 7.40627;  // 2.7 eta
  seedingCfg.sigmaScattering = 50;
  seedingCfg.radLengthPerSeed = 0.1;
  seedingCfg.minPt = 500.;
  seedingCfg.bFieldInZ = 0.00199724;
  seedingCfg.beamPosX = 0;
  seedingCfg.beamPosY = 0;
  seedingCfg.impactMax = 3.;


You can run the seeding example with a comand like this:

.. code-block:: console

   $ <build>/bin/ActsExampleSeedingGeneric \
	--input-dir=data/sim_generic/ttbar_pu200 \
	--output-dir=output_generic_ttbar_pu200 

After running this example, you should see something like this:

.. code-block::
   
   SeedingPerfo  INFO   Efficiency (nMatchedParticles / nAllParticles) = 0.917263
   SeedingPerfo  INFO   Fake rate (nUnMatchedSeeds / nAllSeeds) =0.317763
   SeedingPerfo  INFO   Duplication rate (nDuplicatedMatchedParticles / nMatchedParticles) =0.998711

This output shows the efficiency, fake rate, and duplicate rate for the selected particles.
The example also generates output root files in the output directory.
In ``performance_seeding_hists.root``, you can find the efficiency plots.
The plots below are examples of the efficiency plots produced using the ttbar sample with 200 pile-up vertices.

.. image:: ../figures/performance/fitter/trackeff_vs_eta_ttbar_pu200.png
   :width: 300

.. image:: ../figures/performance/fitter/trackeff_vs_pT_ttbar_pu200.png
   :width: 300


Using a different detector
-------------------
The example above uses the Generic detector, but you can also try the seeding example with a DD4hep detector.
To run the seeding example with the OpenData detector,

.. code-block:: console

   $ <build>/bin/ActsExampleSeedingDD4hep \
    --dd4hep-input ../../acts/Examples/Detectors/DD4hepDetector/compact/OpenDataDetector/OpenDataDetector.xml \
    --input-dir=data/sim_dd4hep/ttbar_pu200 \
    --output-dir=output_dd4_ttbar_pu200

The input dataset needs to be simulated with the same detector in advance.

The detector volumes and layers used for seeding are configured in ``acts/Examples/Run/Seeding/DD4hep/DD4hepSeedingExample.cpp`` as follows.

.. code-block::
		
      // open detector barrel layers
      // the selection intentionally contains duplicates to demonstrate the
      // automatic selection normalization. setting only the volume already
      // selects all layers within it. the explicit layers in the selection
      // should have no effect.
      Acts::GeometryIdentifier().setVolume(13),
      // open detector positive endcap layers
      Acts::GeometryIdentifier().setVolume(14),
      // open detector negative endcap layers
      Acts::GeometryIdentifier().setVolume(12),

If you want to try a different detector geometry, you need to set the detector volumes and layers properly in this file.




