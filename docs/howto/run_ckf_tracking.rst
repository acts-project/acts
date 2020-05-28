Run the CombinatorialKalmanFilter (CKF) tracking example
===============================

Prerequisites
-------------

Acts must be build with activated examples and Pythia8 support
(``ACTS_BUILD_EXAMPLES_PYTHIA8=on``) to enable the fast simulation. ``<build>``
is used to identify the path to the build directory.

Generate a simulation dataset
-----------------------------

Generate a simulation dataset with ttbar process with an average of 200 additional pile-up interactions based on the generic example detector in a 2T magnetic field:

.. code-block:: console

  $ <build>/bin/ActsSimFatrasGeneric \
       --evg-input-type=pythia8 \
       --evg-hard-process Top:qqbar2ttbar=on \
       --evg-pileup=200 \
       --select-pt-gev '0.1:' \
       --select-eta '-2.5:2.5' \
       --fatras-pmin-gev 0.1 \
       --remove-neutral 1 \
       --bf-value=0 0 2 \
       --output-dir=sim_ttbar_pu200 \
       --output-csv=1 \
       --events=10

Setting the output to CSV is necessary since the CKF tracking only reads
CSV files at the moment. 

The above simulation also includes particle selection at different phase, e.g. only generated particles with pT > 100 MeV 
(``--select-pt-gev '0.1:'``) and |eta| <= 2.5 (``--select-eta '-2.5:2.5'``) are passed to simulation.
Further particle selection, e.g. requiring a minimum momentum at 100 MeV (``--fatras-pmin-gev 0.1``) and removing neutral particles (``--remove-neutral 1``), is done during the simulation.

Run the CKF tracking
----------------------

Run the CKF tracking tool that reads the simulation output (truth hits and truth particles), creates smeared
measurements from the true hits, creates seeds (i.e. starting track parameters) from the pre-selected truth particles, 
and run the CKF which will perform the track finding and track fitting simultaneously:

Currently, there are two configurable criteria to select compatible source links on a surface with track parameters in CKF:

* Global maximum chi2 of Kalman filtering. This could be set up via ``--ckf-slselection-chi2max``
* Global maximum number of source links on a surface. This could be set up via ``--ckf-slselection-nmax`` 

.. code-block:: console

   $ <build>/bin/ActsRecCKFTracks \
       --input-dir=sim_ttbar_pu200 \
       --bf-value=0 0 2 \
       --ckf-slselection-chi2max 15 \
       --ckf-slselection-nmax 10 \
       --output-dir=rec_ttbar_pu200

The magnetic field setup should be consistent between simulation and CKF tracking.

Look at the CKF tracking performance
----------------------

The CKF tracking will generate a root file named ``performance_ckf.root`` (the name is currently not configurable via the command line) in the ``output-dir``.
This file includes a few efficiency plots showing the CKF efficiency, fake rate, duplication rate and other plots showing detailed reconstruction info etc.

Example plots to show the CKF efficiency, fake rate and duplication rate for the ttbar sample generated above:

.. image:: ../figures/performance/CKF/trackeff_vs_eta_ttbar_pu200.png
   :width: 300

.. image:: ../figures/performance/CKF/fakerate_vs_eta_ttbar_pu200.png
   :width: 300

.. image:: ../figures/performance/CKF/duplicationrate_vs_eta_ttbar_pu200.png
   :width: 300
