Run the FAst TRAck Simulation
=============================

Prerequisites
-------------

Acts must be build with the ``ACTS_BUILD_EXAMPLES_PYTHIA8=on`` option to
activate building the Fatras simulation. This enables the generic detector. To
be able to run e.g. DD4hep-base detectors, also
``ACTS_BUILD_EXAMPLES_DD4hep=on`` must be set.

In the following, it is assumed that Acts was build as described above in a
separate build directory ``<build>``.

Generate TrackML-like datasets
------------------------------

There exists one simulation executable for each type of detectors. For
TrackML-like datasets, the *Generic* detector must be used. This does not
require any additional settings.

The fast simulations generates the initial particles of the truth event and then
propagates them through the detector. At each surface, particle interactions are simulated using parametric models and potential hits are registered.

There are two options for generating the initial particles: a particle gun or
the Pythia8 generator. The following command generates 100 events using the
particle gun, in a 2T magnetic field, and writes the generated events into the
``sim-pg`` directory:

.. code-block:: console

   $ <build>/bin/ActsSimFatrasGeneric \
       --output-dir=sim-gun \
       --output-csv=1 \
       --evg-input-type=gun \
       --bf-value=0 0 2 \
       --events=100

For each event, the following files will be created

-   ``event<number>-particles_generated.csv``
-   ``event<number>-particles_initial.csv``
-   ``event<number>-particles_final.csv``
-   ``event<number>-truth.csv``
-   ``event<number>-hits.csv``
-   ``event<number>-cells.csv``

where ``<number>`` is the event number. The first three files contain the
generated particle states, initital and final states of simulated particles. The
simulated particles differ from the generated particles: particles could not be
simulated due to kinematic cuts or additional particles might have been
generated due to interactions. The truth contains the true intersection with all
surfaces, while the hits and the cells describe the simulated detector readout.

When generating particles with Pythia8, the event is build from a hard scatter
interaction and additional pileup interactions. The following commands will
generate top-pairs as the hard scatter interaction and uses the default soft QCD
process for the additional 140 (on average) pileup interactions. The detector
setup is unchanged.

.. code-block:: console

   $ <build>/bin/ActsSimFatrasGeneric \
       --output-dir=sim-pythia8 \
       --output-csv=1 \
       --evg-input-type=pythia8 \
       --evg-hard-process=Top:qqbar2ttbar=on \
       --evg-pileup=140 \
       --bf-value=0 0 2 \
       --events=100

The output file structure will be the same as discussed above.