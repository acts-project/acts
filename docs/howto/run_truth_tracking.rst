Run the truth tracking examples
===============================

Prerequisites
-------------

Acts must be build with activated examples and Pythia8 support
(``ACTS_BUILD_EXAMPLES_PYTHIA8=on``) to enable the fast simulation. ``<build>``
is used to identify the path to the build directory.

Generate simulation datasets
-----------------------------

Generate two different example simulation datasets based on the generic example detector in a 2T magnetic field:

-  10000 (``--pg-nparticles 10000``) single muons with pT at [0.1, 100 ) GeV (``--pg-pt-range 0.1 100``) and |eta| <= 2.5 (``--pg-eta-range -2.5 2.5``)

.. code-block:: console

   $ <build>/bin/ActsSimFatrasGeneric \
       --evg-input-type=gun \
       --pg-nparticles 10000 \
       --pg-pt-range 0.1 100 \
       --pg-eta-range -2.5 2.5 \
       --bf-value=0 0 2 \
       --output-dir=sim_singlemuon \
       --output-csv=1 \
       --events=10

-  ttbar process with an average of 200 additional pile-up interactions (``--evg-pileup=200``)
The particles are selected at different phase, e.g. only generated particles with pT > 100 MeV 
(``--select-pt-gev '0.1:'``) and |eta| <= 2.5 (``--select-eta '-2.5:2.5'``) are passed to simulation.
Further particle selection, e.g. removing neutral particles (``--remove-neutral 1``), is done during the simulation.

.. code-block:: console

  $ <build>/bin/ActsSimFatrasGeneric \
       --evg-input-type=pythia8 \
       --evg-hard-process Top:qqbar2ttbar=on \
       --evg-pileup=200 \
       --select-pt-gev '0.1:' \
       --select-eta '-2.5:2.5' \
       --remove-neutral 1 \
       --bf-value=0 0 2 \
       --output-dir=sim_ttbar_pu200 \
       --output-csv=1 \
       --events=10

Setting the output to CSV is necessary since the truth tracking only reads
CSV files at the moment. 

Run the truth tracking
----------------------

Run the truth tracking tool that reads the simulation output (truth hits and truth particles), creates smeared
measurements from the true hits, creates seeds (i.e. starting track parameters) from the truth particles, builds truth tracks (i.e. uses the truth
information to group simulated hits into tracks) and fits them. Examples of truth tracking with the different simulation output generated above:

-   Single muon sample

.. code-block:: console

   $ <build>/bin/ActsRecTruthTracks \
       --input-dir=sim_singlemuon \
       --bf-value=0 0 2 \
       --output-dir=rec_singlemuon

-  ttbar sample

.. code-block:: console

   $ <build>/bin/ActsRecTruthTracks \
       --input-dir=sim_ttbar_pu200 \
       --bf-value=0 0 2 \
       --output-dir=rec_ttbar_pu200

The magnetic field setup should be consistent between simulation and truth tracking. 

Look at the truth tracking performance
----------------------

The truth tracking will generate three root files (the name of those root files are currently not configurable via the command line) in the ``output-dir``:

*   ``tracks.root``
This includes a tree with one entry representing one trajectory. From this file, one could check the information of every measurement track state on the trajectory.

*  ``performace_track_finder.root``
This includes a tree showing performance of the truth track finding.

*  ``performance_track_fitter.root``
This includes a few histograms showing the residual and pull of the fitted perigee track parameters and efficiency plots showing the fitting efficiency etc.

Example plots to show the fitting efficiency versus eta and pT for ttbar sample generated above:

.. image:: ../figures/performance/fitter/trackeff_vs_eta_ttbar_pu200.png
   :width: 300

.. image:: ../figures/performance/fitter/trackeff_vs_pT_ttbar_pu200.png
   :width: 300

Example plots to show the average number of measurments and holes versus eta for ttbar sample generated above:

.. image:: ../figures/performance/fitter/nMeasurements_vs_eta_ttbar_pu200.png
   :width: 300

.. image:: ../figures/performance/fitter/nHoles_vs_eta_ttbar_pu200.png
   :width: 300

To draw the resolution (residual and pull) of fitted perigee track parameters for e.g. ttbar sample, one could use:

.. code-block:: console

 $ root <source>/Examples/Scripts/perigeeParamResolution.C("rec_ttbar_pu200/performance_track_fitter.root")'

``<source>`` here is used to identify the path of the source directory. 

An example plot of the pull distribution of fitted perigee track parameters for the ttbar sample generated above:

.. image:: ../figures/performance/fitter/pull_perigee_parameters_ttbar_pu200.png
   :width: 600
