.. _simulate-TrackML:

Run the FAst TRAck Simulation
=============================

Prerequisites
-------------

Acts must be build with the ``ACTS_BUILD_EXAMPLES=on`` option to produce the
Fatras simulation executables. With this option alone, Fatras can be run for
the generic detector using a particle gun or externally generated particle
input data.

Additional options might be needed to enable different generators or detectors.
The ``ACTS_BUILD_EXAMPLES_PYTHIA8=on`` option enables Pythia8-based event
generator executables. To be able to run e.g. DD4hep-base detectors, the
``ACTS_BUILD_EXAMPLES_DD4hep=on`` option must be set. The full list of
available options is available in the *Getting started* guide.

In the following, it is assumed that Acts was build with Pythia8 support as
described above in a separate build directory ``<build>``.

Generate particle truth information
-----------------------------------

Before running the fast simulation, a truth dataset needs to be generated. The
simplest option is to use the particle gun, which generates particles
uniformly distributed within the given ranges. The particle gun does not
implement a physics process but is useful e.g. single particle studies.

The following command will generate 100 events with four (anti-)muons in
each one and write the data as CSV files into the ``data/gen/four_muons``
directory.

.. code-block:: console

   $ <build>/bin/ActsExampleParticleGun \
       --events=100 \
       --output-dir=data/gen/four_muons \
       --output-csv \
       --gen-phi-degree=0:90 \
       --gen-eta=-2:2 \
       --gen-mom-gev=1:5 \
       --gen-pdg=13 \
       --gen-randomize-charge \
       --gen-nparticles=4

The muons are generated within the given kinematic range, but are independent
from each other. With the particle gun, all particles always originate from
a single common vertex. With the default settings, the vertex will always be
at the origin. Additional options allow to smear the vertex positions. The
full list of options is available via

.. code-block:: console

   $ <build>/bin/ActsExampleParticleGun -h

To generate more realistic truth datasets you can use the Pythia8-based
generator. When generating particles with Pythia8, the event is build from
a fixed number of hard scatter interactions and additional pileup interactions.
The number of pileup interactions is drawn from a Poisson distribution.
The following commands will generate proton-proton collisions with a
center-of-mass energy of 14TeV that produce a top-pair from the single hard
scatter interaction and additional 140 (on average) pileup interactions using
the default soft QCD process. The output is written as CSV files in the
``data/gen/ttbar_mu140`` directory.

.. code-block:: console

   $ <build>/bin/ActsExamplePythia8 \
       --events=100 \
       --output-dir=data/gen/ttbar_mu140 \
       --output-csv \
       --rnd-seed=42 \
       --gen-cms-energy-gev=14000 \
       --gen-hard-process=Top:qqbar2ttbar=on \
       --gen-npileup=140

A full list of options and their default values is available via

.. code-block:: console

   $ <build>/bin/ActsExamplePythia8 -h

Simulate a TrackML-like detector
--------------------------------

There is one simulation executable for each type of detectors. For
TrackML-like datasets, the *Generic* detector must be used. The simulation
requires particle truth information as input. This can be either generated
on the fly using a particle gun or read from files that were created
separately (as discussed above). For each input particle, the fast
simulations propagates it through the detector. At each surface, particle
interactions are simulated using parametric models and potential hits are
registered.

The following command simulates 100 events for the generic detector in a 2T
magnetic field and writes the generated events as CSV files into the
``data/sim_generic/single_muon`` directory. The particle truth is
generated on-the-fly using the integrated particle gun; by default single
muons in a reasonable kinematic range are generated.

.. code-block:: console

   $ <build>/bin/ActsExampleFatrasGeneric \
       --output-dir=data/sim_generic/single_muon \
       --output-csv \
       --events=100 \
       --bf-constant-tesla=0:0:2

For each event, the following files will be created

-   ``event<number>-particles_initial.csv``
-   ``event<number>-particles_final.csv``
-   ``event<number>-hits.csv``

where ``<number>`` is the event number. The first two files contain the
initital and final states of simulated particles. The
simulated particles can differ from the generated input particles: particles might not have been
simulated due to kinematic cuts or additional particles might have been
generated due to interactions. The hits describe the simulated detector readout.

To use some of the previously generated truth datasets, the ``--input-dir``
option must be set. The following command reads the previously generated
top-pair sample with some additional selections cuts on the generated
particles.

.. code-block:: console

   $ <build>/bin/ActsExampleFatrasGeneric \
       --input-dir=data/gen/ttbar_mu140 \
       --output-dir=data/sim_generic/ttbar_mu140 \
       --output-csv \
       --select-eta=-3:3 \
       --select-pt=0.5: \
       --remove-neutral \
       --bf-constant-tesla=0:0:2 

The output file structure will be the same as discussed above.

Simulate a DD4hep based detector
--------------------------------

Similar to the simulation of the generic detector, another example detector called *Open Data Detector* (ODD) based on the DD4hep description can be simulated with the ``ACTS_BUILD_DD4HEP_PLUGIN`` option enabled. For the DD4hep based detector, a xml file for the geometry description must be specified by the option ``--dd4hep-input``. In addition, a material mapping file must be provided via the option ``--mat-input-file`` in order to take into account the material effects in simulation. 

The following command simulates 100 single muon events for the ODD in a 2T magnetic field.

.. code-block:: console

   $ <build>/bin/ActsExampleFatrasDD4hep \
       --dd4hep-input <source>/thirdparty/OpenDataDetector/xml/OpenDataDetector.xml \ 
       --mat-input-type file \
       --mat-input-file <source>/thirdparty/OpenDataDetector/config/odd-material-mapping.config \
       --output-dir=data/sim_ODD/single_muon \
       --output-csv \
       --events=100 \
       --bf-constant-tesla=0:0:2

The following command simulates the top-pair events based on the previously generated
top-pair sample for the ODD in a 2T magnetic field.

.. code-block:: console

   $ <build>/bin/ActsExampleFatrasDD4hep \
       --dd4hep-input <source>/thirdparty/OpenDataDetector/xml/OpenDataDetector.xml \ 
       --mat-input-type file \
       --mat-input-file <source>/thirdparty/OpenDataDetector/config/odd-material-mapping.config \
       --input-dir=data/gen/ttbar_mu140 \
       --output-dir=data/sim_generic/ttbar_mu140 \
       --output-csv \
       --select-eta=-3:3 \
       --select-pt=0.5: \
       --remove-neutral \
       --bf-constant-tesla=0:0:2
