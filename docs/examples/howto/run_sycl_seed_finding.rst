Run the SYCL seed finding tests
===============================

List available devices
----------------------

We can list available platforms and devices with **clinfo**:

.. code-block:: bash

    [bash][atspot01]:build > clinfo
        Number of platforms                               3
        Platform Name                                   Intel(R) OpenCL
        Platform Vendor                                 Intel(R) Corporation
        Platform Version                                OpenCL 2.1 LINUX
        ...

This is useful to decide which device we want to generate binaries for (because we need to give that to CMake). Note, that this would not only list OpenCL, but also CUDA devices.

We can install *clinfo* from our packet manager.

Compiler
--------

The SYCL seed finding implementation includes extensions from DPC++, so the code should be compiled with a dpcpp compiler. It can be cloned and built from the Intel LLVM `github`_ repository, or it could also be downloaded and installed from the `oneAPI`_ website.

.. _github: https://github.com/intel/llvm/
.. _oneAPI: https://software.intel.com/content/www/us/en/develop/articles/intel-oneapi-dpcpp-compiler-release-notes-beta.html

CMake should inform us whether our compiler is SYCL compatible.

.. code-block::

    -- Checking if /atlas/software/intel/clang/12.0.0-2020-08-24/x86_64-ubuntu1804-gcc8-opt/bin/clang++ is SYCL capable... success
    -- Checking for available SYCL target(s)...
    --   - Found target: nvptx64-nvidia-cuda-sycldevice
    -- Checking for available SYCL target(s)... done
    ...

Building with CMake
-------------------

During the CMake configuration the flag ``ACTS_BUILD_PLUGIN_SYCL=on`` needs to be set. CMake will check whether our compiler is SYCL compatible and will look for available SYCL target architectures.

By default, SPIR 64-bits targets are looked for, but this can be configured with the ``SYCL_POSSIBLE_TARGETS`` variable.
For example, to target CUDA backends, we should set this variable to ``nvptx64-nvidia-cuda-sycldevice``.

.. code-block:: cmake

    # Figure out which SYCL target platforms are available.
    set( SYCL_POSSIBLE_TARGETS "spir64-unknown-unknown-sycldevice"
    ...

Also, by setting the variable ``SYCL_SEARCH_SUPPORT_LIBRARIES`` it is possible to specify support libraries or object files that are needed for execution on a specific target. Missing object files cause a runtime error.

Trying the code
---------------

We can build the unit test ActsUnitTestSeedFinderSycl for this purpose (actually, this is not a unit test).
The test compares the results of the CPU and SYCL seed finding algorithm in terms of speed and precision.

It takes the following command-line options:

.. code-block:: bash

    Allowed options:
    -h [ --help ]           Print usage message.
    -f [ --FILE ] arg       Provide path for input file.
    -n [ --NUM ] arg (=500) Number of groups to iterate in seed finding.
    -d [ --DEVICE ] arg     Provide a substring of the preferred device.
    -l [ --LIST ]           List available SYCL platforms and devices.
    -G [ --GPU ]            Execute code only on gpu. Default is 0.
    -a [ --ALL ]            Analyze all groups. Default is 0.
    -m [ --MATCH ]          Count seed matches. Default is 0.

We can select the preferred device by the -d option, or we could also set the ``SYCL_BE`` environment variable.

By setting the environment variable ``SYCL_BE`` we can alter the behavior of the SYCL ``default_selector``
For example by setting it to ``PI_CUDA`` it forces the usage of the CUDA backend (if available).
Similarly ``PI_OPENCL`` forces the usage of the OpenCL backend.

We have to provide an input file with the -f flag.
The input file should have a predefined format:

.. code-block:: bash

    [bash][atspot01]:build > head /atlas/acts_data/atlas_seeds/pu200/evt10.txt 
        lxyz 1 -35.4702 -10.973 -346 0.05 0.1
        lxyz 1 -37.5198 -12.2698 -383 0.05 0.1
        lxyz 1 -39.7144 -13.9292 -423 0.05 0.1
        lxyz 1 -41.723 -16.0801 -470 0.05 0.1
        lxyz 1 -44.332 -18.1633 -525 0.05 0.1

Numbers correspond to layer, x, y, z coordinates, R, and Z variances (?).

In the end, we should see something like this:

.. code-block:: bash

    [bash][atspot01]:build > SYCL_BE=PI_CUDA bin/ActsUnitTestSeedFinderSycl -f /atlas/acts_data/atlas_seeds/pu200/evt10.txt -m
        11:48:57    QueueWrapper   INFO      Running on: GeForce RTX 2060
        read 360734 SP from file /atlas/acts_data/atlas_seeds/pu200/evt10.txt
        Preparation time: 0.637518
        Analyzed 260 groups for CPU
        Analyzed 260 groups for SYCL

        ------------------------- Time Metric -------------------------
                     Device:        CPU       SYCL  Speedup/ Agreement
                   Time (s):  13.316532   1.733498            7.681883
                Seeds found:     171516     171516           99.950447
        ---------------------------------------------------------------
