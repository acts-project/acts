ACTS Material Mapping Auto-Tuning Tutorial
===========================================

(This how to uses **deprecated** executables and will be updated soon)

The material mapping can be quite a cumbersome process, especially when used on a detector with an evolving geometry. The point of the auto-tuning is to replace the person power needed to perform the optimisation by computing power and time.  This tutorial will present how to perform the auto-tuning of the material mapping process. It assumes the reader is familiar with the concepts related to the material mapping. For more information, refer yourself to the `run_material_mapping` how to guide. This page will explain how to perform the auto-tuning material mapping with the Open Data Detector (ODD) but it can be applied to any other detector (see the other how to change detector).

Prerequisites
--------------
The prerequisites are the same as for the material mapping, we will need the ACTS Examples, Geant4 and the JSON plugin (``ACTS_BUILD_EXAMPLES``, ``ACTS_BUILD_EXAMPLES_GEANT4`` and ``ACTS_BUILD_PLUGIN_JSON``). Depending on the type of detector we want to map, you will need to use some additional packages, in our case ``ACTS_BUILD_EXAMPLES_DD4HEP`` and ``ACTS_BUILD_PLUGIN_TGEO`` are needed. The auto-tuning is implemented as part of the Acts python bindings, we will thus also need to have them installed (``ACTS_BUILD_EXAMPLES_PYTHON_BINDINGS``).

We will also need an optimisation library. The material mapping auto-tuning has been implemented using Orion (https://orion.readthedocs.io/en/stable/), it can be easily installed using the following pip command:

.. code-block:: console

   $ pip install orion

Configuration
--------------

Similarly to the regular material mapping, we first need to configure which surfaces will be optimised in the mapping. For this we first need to extract the list of all the surfaces in our detector, to do so we will use the GeometryExample:

.. code-block:: console

   $ <build>/bin/ActsExampleGeometryDD4hep -n1 -j1 \
       --mat-output-file geometry-map \
       --dd4hep-input <source>/thirdparty/OpenDataDetector/xml/OpenDataDetector.xml \
       --output-json \
       --mat-output-allmaterial true \
       --mat-output-sensitives false

The JSON file can then be edited to select which surfaces we want to have material mapped on. To add one surface to the material mapping, one simply needs to switch the ``mapMaterial`` variable to ``true``. The binning doesn't need to be changed, as it will be decided by the auto-tuning algorithm itself.

As a rule of thumb, the material can be mapped on the approach surface of your sensitive surface and the boundary between the different volume of your detector.

.. warning::
  When performing the auto-tuning, the list of surface the material is mapped onto must stay consistent. You want to change the surface used, you will need to restart the tuning from the start. 

Geantino scan
--------------

The next step is to do a geantino scan of our detector. For this, we will use the ``MaterialRecording`` application:

.. code-block:: console

   $ <build>/bin/ActsExampleMaterialRecordingDD4hep -n1000 -j1 \
       --dd4hep-input <source>/thirdparty/OpenDataDetector/xml/OpenDataDetector.xml \
       --output-root


The result of the geantino scan will be a root file containing material tracks. Those contain the direction and production vertex of the geantino, the total material accumulated and all the interaction points in the detector.

Material Mapping Optimisation 
------------------------------

Once we have our json configuration file and our geantino material track file we can start the optimisation of the material map. For this, we will use the ``material_mapping_optimisation.py`` python binding file. This script file has five options that can be changed:

- ``numberOfJobs`` the number of simultaneous jobs executed (this can be as high as you want, but if you have enough Core 40 works quite well).
- ``topNumberOfEvents`` the number of events from the input material track file used in the mapping can be changed (by default: 10000) it needs to be the same in all the calls to ``material_mapping_optimisation``.
- ``inputPath`` the path to the input directory
- ``outputPath`` the path to the output directory
- ``doPloting`` if added at the end of the optimisation script, the optimal material will be computed and result plots will be created. (Should be used at the end of the optimisation).
- ``readCachedSurfaceInformation`` if added the material-surface association will be taken from the input material track file (doesn't work with geantino file, you need to use the material track file obtained from running the material mapping).

When using Orion our optimisation algorithm can easily be parallelised. Orion uses a Database system that stores all the mapping configuration tested so far, before starting a new trial the database is checked that the parameters have not been used before. Each time you run ``material_mapping_optimisation.py`` a batch of ``numberOfJobs`` trials will be performed and added to the database. Once you have run it enough time, you can extract the result by running with the ``--doPloting`` option.

The current python script uses a simple random search, more performant algorithms could be tried. You can refer to the Orion documentation to see how to implement/uses them.

To simplify the use of the auto-tuning a bash script is available in ``Examples/Scripts/Python/Auto-tuning/Orion/``, performing the optimisation should be as simple as launching this script. You will only need to modify the path to your input and output directory that are called ``MaterialMappingInputDir`` and ``MaterialMappingOutputDir`` by default.

Material Validation
--------------------

Once you have run the optimisation enough time and extracted the optimised material map, you can run material map validation from the ``run_material_mapping`` to validate your map. 


Implementation 
---------------

This section will present in more detail how the optimisation is implemented in ``material_mapping_optimisation.py``. This information is not necessarily needed to run it, but might be useful if you want to modify the script.

When calling ``material_mapping_optimisation.py``, we create one process per surface using python ``multiprocessing``. Each of those trial processes will be in charge of optimising the binning for the corresponding surface and connect to the corresponding database. This is needed because Orion can only have one active database per process. In each of those processes we will create ``numberOfJobs`` trials, each trial corresponding to a different binning, then pipe those binning to the main process.

After creating the trials process wait until it has received one binning per trial process (so one per surface), those are combined to configure one material mapping job that will be launched in a separate process. This is performed a number of times equal to the value of ``numberOfJobs``.

Each of the mapping process will need to perform the material mapping twice, once to determine the average material in each bin and a second time to compute the variance. For each surface, a score is then computed using the variance and the number of hits in each bin. This score is then piped to the main process.

After receiving the scores, the main process pipes them back to the trial processes which will store them in the database. The script then ends when this has been done for all the jobs. If the ``--doPloting`` option was used, each trial process will also return some plot related to the optimisation performance so far. They will also pipe their best binning to the main process. One last mapping job is then performed in the main process, this will return the optimised material map and a material track file where the material is already associated to the surfaces. This last file can be used for validation and as an input to futur mapping jobs to speed them by up to 50% (using the ``--readCachedSurfaceInformation`` option).

![Diagramme of the material mapping auto-tuning](/figures/materialMapping/ActsMaterialMappingAutoTuning.png)