ACTS Tutorial on Auto-Tuning in CombinatorialKalmanFilter (CKF)
===============================================================
The tracking algorithms require a number of pre-initialized parameters that are often hand-tuned to obtain high performance. Usually, the value of these parameters change as the underlying geometrical or magnetic configuration changes. An automatic tuning of these parameters can be very useful for obtaining highly efficient parameter configuration as well as for studying different detector geometries. This tutorial is based on parameter optimization studies using two different optimization frameworks: Optuna and Orion. Eight parameters of Track Seeding algorithm have been tuned using these frameworks and their performance have been studied on CKF.

Prerequisites
-------------
Since Optuna and Orion are independent frameworks, these need to be installed separately. The following setup works on any machine with cvmfs access that requires creating a virtual python environment:

.. code-block:: console

   $ source /cvmfs/sft.cern.ch/lcg/views/LCG_100/x86_64-centos7-gcc10-opt/setup.sh
   $ python3 -m venv PYTHON_VIRTUAL_ENV
   $ source PYTHON_VIRTUAL_ENV/bin/activate
   $ export PYTHONPATH=
   $ python -m pip install --upgrade pip
   $ pip install -r acts/Examples/Python/tests/requirements.txt
   $ pip install pytest --upgrade
   $ pip install cmake dataclasses sphinxcontrib-applehelp sphinxcontrib-jsmath sphinxcontrib-serializinghtml argparse sphinxcontrib-devhelp sphinxcontrib-htmlhelp sphinxcontrib-qthelp AppDirs filelock joblib pandas plotly psutil pyYAML requests scipy tabulate cloudpickle scikit-learn orion==0.2.2 cloudpickle==1.6.0 optuna matplotlib

Once all the dependencies are installed, Build ACTS with python bindings and Pythia8 support.

.. code-block:: console

   $ cmake -DACTS_BUILD_EXAMPLES_PYTHIA8=ON -DACTS_BUILD_PLUGIN_DD4HEP=OFF -DACTS_BUILD_PLUGIN_JSON=ON -DACTS_BUILD_PLUGIN_ROOT=ON -DACTS_BUILD_EXAMPLES_DD4HEP=OFF -DACTS_BUILD_EXAMPLES_GEANT4=ON -DACTS_BUILD_INTEGRATIONTESTS=OFF -DACTS_BUILD_UNITTESTS=OFF -DACTS_BUILD_EXAMPLES_PYTHON_BINDINGS=ON -DACTS_BUILD_ODD=OFF -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_STANDARD=17 -S . -B build/
   $ make
   $ source build/python/setup.sh

Once this setup is ready, at each new login, just do:

.. code-block:: console

   $ source /cvmfs/sft.cern.ch/lcg/views/LCG_100/x86_64-centos7-gcc10-opt/setup.sh
   $ source PYTHON_VIRTUAL_ENV/bin/activate
   $ export PYTHONPATH=
   $ source build/python/setup.sh

How auto-tuning works
---------------------
A list of parameters and their range are provided to optimization framework. The framework predict a value for each parameter and provide them to the tracking algorithm. The tracking algorithm runs with these values and return the output to the optimization framework. The optimization framework analyzes the output and provide a new set of values to tracking algorithm and this cycle continues for a number of trials. At each trial, the optimization framework computes a score based on the tracking algorithm output and tries to minimize/maximize this score by suggesting better parameter configurations. The current frameworks compute the score based on the CKF track efficiency, fake rate, duplicate rate and run-time:

Score = Efficiency - (fakeRate + DuplicateRate/k_dup + run-time/k_time)

where k_dup and k_time are the weights for duplicate rate and run-time. These weights play a significant role in determining the configuration of best performing parameters from the optimization frameworks.

The list of track seeding parameters that are auto-tuned using these frameworks are as follows:
* maxPtScattering: upper p_T limit for scattering angle calculations
* impactMax: maximum value of impact parameter
* deltaRMin: minimum distance in r between two measurements within one seed
* deltaRMax: maximum distance in r between two measurements within one seed
* sigmaScattering: number of sigma used for scattering angle calculations
* radLengthPerSeed: average radiation lengths of material on the length of a seed
* maxSeedsPerSpM: number of space-points in top and bottom layers considered for compatibility with middle space-point
* cotThetaMax: maximum cotTheta between two space-points in a seed to be considered compatible

The python scripts utilizing these optimization using Optuna and Orion frameworks are present in acts/Examples/Scripts/Optimization directory.

Run auto-tuning using Optuna
----------------------------
The Optuna auto-tuning script can be run directly by invoking:
`` python Optuna_tuning.py``

This creates a new optuna study for a given number of trials defined within the script. The direction is set to maximize which means that the framework will try to maximize the score.

.. code-block:: console

   $ study = optuna.create_study(study_name=study_name,
		storage="sqlite:///{}.db".format(study_name),
		direction='maximize',
		load_if_exists=True)

The objective function defines the list of parameters to tune along with their range. It suggests a parameter configuration using ``trial.suggest_float`` or ``trial.suggest_int`` function of optuna and runs CKF for this parameter configuration. The output of CKF is read from the performance file named ``performance_ckf.root`` and a score function is constructed within the objective function.

The objective function and the number of trials are passed to optuna optimize function:
``study.optimize(objective, n_trials=100)``. The best parameter configuration after all the trials is read from ``study.best_trial.params.items()``.

Run auto-tuning using Orion
---------------------------
The Orion auto-tuning script can be run directly by invoking:
`` python Orion_tuning.py``

A dictionary called space is created by providing a list of parameters and their range and new experiment is build using

.. code-block:: console

   $ experiment = build_experiment(
		"orion_new",
		space=space,
		storage=storage)

The objective function picks up a value for each parameter, run CKF and construct a score function from CKF output as in Optuna case. The only difference is that it tries to minimize the score unlike optuna. The objective function and number of trials are passed to the orion workon function
``experiment.workon(objective, max_trials=100)``.

The best parameter configuration corresponds to the minimum score function value and can be obtained from the experiment.
