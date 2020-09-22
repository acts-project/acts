SYCL plugin
===========

The SYCL plugin allows building interfaces that implement algorithms in SYCL.

Build
-----

During the CMake configuration the flag ``ACTS_BUILD_PLUGIN_SYCL=on`` needs to be set. CMake will check whether our compiler is SYCL compatible and will look for available SYCL target architectures.
By default, SPIR 64-bits targets are looked for, but this can be configured with the ``SYCL_POSSIBLE_TARGETS`` variable.
For example, to target CUDA backends, we should set this variable to ``nvptx64-nvidia-cuda-sycldevice``.

Also, by setting the variable ``SYCL_SEARCH_SUPPORT_LIBRARIES`` it is possible to specify support libraries or object files that are needed for execution on a specific target. Missing object files cause a runtime error.

Device selection
----------------

Firstly, we can list available platforms and devices with **clinfo**:

.. code-block:: bash

    [bash][atspot01]:build > clinfo
        Number of platforms                             3
        Platform Name                                   Intel(R) OpenCL
        Platform Vendor                                 Intel(R) Corporation
        Platform Version                                OpenCL 2.1 LINUX
        ...

CMake should recognize these platforms during configuration (if not told otherwise).

There are more options we can choose from to configure device selection (meaning constructing a ``cl::sycl::device_selector`` object).

We can alter the behavior of the SYCL ``default_selector`` by setting the environment variable ``SYCL_BE`` (BE stands for backend).
For example by setting it to ``PI_CUDA`` it forces the usage of the CUDA backend (if available).
Similarly ``PI_OPENCL`` forces the usage of the OpenCL backend. For more information see this Getting Started `guide`_.:

.. code-block:: bash

    SYCL_BE=PI_CUDA <binary> <arguments>

We can also use a custom device selector. As an example, see :class:`Acts::Sycl::DeviceSelector`, which selects a non-OpenCL CUDA backend.
This is the device selector used by the :class:`Acts::Sycl::QueueWrapper` object that manages the construction, ownership and destruction of a ``cl::sycl::queue`` object.
This solution allows us to construct a ``queue`` object in translation units that are not linked against SYCL.

.. _guide: https://intel.github.io/llvm-docs/GetStartedGuide.html#run-simple-dpc-application

Algorithms
----------

Currently, the following algorithms are implemented in SYCL:

  #. Seeding

SYCL Seeding
------------

Usage
^^^^^

To create a Seedfinder object that implements the seeding algorithm in SYCL, we need to instantiate an object from :class:`Acts::Sycl::Seedfinder`.
Similarly to :class:`Acts::Seedfinder`, we need to provide a configuration object of :class:`Acts::SeedfinderConfig<external_spacepoint_t>` and an object that implements *experiment specific* cuts.
As these cuts are performed in SYCL kernels (on the device side), instead of a :class:`IExperimentCuts` instance, we need to construct a :class:`DeviceExperimentCuts` one.

.. code-block:: cpp

    Seedfinder(
      Acts::SeedfinderConfig<external_spacepoint_t> config,
      const Acts::Sycl::DeviceExperimentCuts& cuts,
      Acts::Sycl::QueueWrapper wrappedQueue = Acts::Sycl::QueueWrapper());

In the current implementation, the member functions :func:`DeviceExperimentCuts::seedWeight()` and :func:`DeviceExperimentCuts::singleSeedCut` in the header file ``DeviceExperimentCuts.hpp`` need to be rewritten to have our custom experiment cuts.

.. code-block:: cpp

    float seedWeight(const detail::DeviceSpacePoint& bottom,
                     const detail::DeviceSpacePoint& middle,
                     const detail::DeviceSpacePoint& top) const {...}
    /*...*/

    bool singleSeedCut(float weight, const detail::DeviceSpacePoint& bottom,
                        const detail::DeviceSpacePoint& middle,
                        const detail::DeviceSpacePoint& top) const {...}

Optionally we can also give a :class:`Acts::Sycl::QueueWrapper` object to the constructor of :class:`Acts::Sycl::Seedfinder`, which is a wrapper object around a ``cl::sycl::queue`` type.
This allows us to construct our own ``queue`` instance and to reuse it.

Implementation details
^^^^^^^^^^^^^^^^^^^^^^

The following section describes memory management, kernel scheduling and array indexing for the SYCL seed finding algorithm.

We start out with the duplet search that looks for compatible bottom-middle and middle-top space point pairs.

In case we have 5 middle SP and 4 bottom SP, our temporary array of
the compatible bottom duplet indices would look like this:

.. figure:: ../figures/plugins/sycl/duplet_search_matrix.png
  :name: duplet_search_matrix
  :align: center
  :width: 200

  Rows correspond to middle space points, numbers are bottom space point indices. Threads are executed concurrently, so the order of bottom SP indices is random.

We will flatten this matrix out, and store the indices the following way:

.. figure:: ../figures/plugins/sycl/flat_matrix.png
  :name: flat_matrix
  :align: center
  :width: 350

  Storing bottom SP indices for all middle SPs.

To be able to get the indices of middle SPs in constant time inside kernels, we will also prepare arrays that store the indices of the middleSPs of the edges.

(For the same purpose, we could also do a binary search on the array on :numref:`prefix_sum_array`, and we will do exactly that later, in the triplet search kernel.)

.. figure:: ../figures/plugins/sycl/middle_duplet_indices.png
  :name: middle_duplet_indices
  :align: center
  :width: 350

To find out where the indices of bottom SPs start for a particular middle SP, we use prefix sum arrays.
We know how many duplets were found for each middle SP:

.. figure:: ../figures/plugins/sycl/count_duplets.png
  :name: count_duplets
  :align: center
  :width: 250

We will make a prefix sum array of these counts, with a leading zero:

.. figure:: ../figures/plugins/sycl/prefix_sum_array.png
  :name: prefix_sum_array
  :align: center
  :width: 300

  Prefix sum array of the counted values of compatible bottom SPs per middle SP.

If we have the middle SP with index 1, then we know that the indices of the compatible bottom SPs are in the range (left closed, right open) [2,5) of the previously flattened array in :numref:`flat_matrix`.
In this case, these indices are 3 and 2, so we'd use these to index deviceBottomSPs to gather data about the bottom SP.


In this example, will execute the coordinate transformation on 7 threads.

The size of the array storing our transformed coordinates is also 7, the sum of bottom duplets we found so far.

The process for middle-top space points is the same.

For the triplet search, we calculate the upper limit of constructible triplets.

For this, we multiply the number of compatible bottom and compatible top SPs for each middle SP, and add these together. This is 
:math:`nb_0*nt_0 + nb_1*nt_1 + ...` where :math:`nb_k` is the number of compatible bottom SPs for the :math:`k`th middle SP, similarly :math:`nt_k` is for tops.

We construct a prefix sum array (of length :math:`M+1`) of the calculated combinations. (Where :math:`M` is the number of middle space points.)

.. _table-1:

.. list-table:: Max triplet combinations prefix sum array
   :widths: 25 25 25 10 25
   :header-rows: 1

   * - middle SPs
     - :math:`middleSP_0`
     - :math:`middleSP_1`
     - ...
     - :math:`middleSP_M`

   * - number of combinations
     - :math:`nb_0*nt_0`
     - :math:`nb_0*nt_0 + nb_1*nt_1`
     - ...
     - :math:`\sum_{i=0}^{M}nb_i*nt_i`

We will start kernels and reserve memory for these combinations but only so much we can fit into memory at once.

For later, let :math:`MAM` be the maximum allocatable memory for triplet search.

We start by adding up summing the combinations, until we arrive at a :math:`k` which for:

.. math::

    \sum_{i=0}^{k+1}nb_i*nt_i > MAM

(or :math:`k = M`).

So we know, that we need to start our first kernel for the first :math:`k` middle SPs.

Inside the triplet search kernel we start with a binary search, to find out which middle SP the thread corresponds to.
Note, that the array in `table-1`_ is a monotone increasing series of values which allows us to do a binary search on it.

Inside the triplet search kernel we count the triplets for fixed bottom and middle SP.

The triplet filter kernel is calculated on threads equal to all possible bottom-middle combinations for the first :math:`k` middle SPs, which are the sum of bottom-middle duplets.

If the triplet search and triplet filter kernel finished, we continue summing up possible triplet combinations from the :math:`(k+1)` th middle SP.

Inside the kernels we need to use an offset, to be able to map threads to space point indices.

Performance
^^^^^^^^^^^

It depends on our architecture, the size of the event we are reconstructing, and the effectiveness of our experiment specific cuts how well the algorithm performs, and whether we can benefit at all from using the SYCL plugin.
It is advised to compare performance and precision first with the CPU version of the chosen algorithm.
This should be possible with the tests provided.

Resources
---------

For more information about SYCL see the `specification`_ (date: 2020. September 7.).
There is a `documentation`_ for Intel implementation and DPC++ extensions (see `examples`_).

.. _documentation: https://software.intel.com/content/www/us/en/develop/download/intel-oneapi-programming-guide.html
.. _specification: https://www.khronos.org/registry/SYCL/specs/sycl-1.2.1.pdf
.. _examples: https://github.com/intel/llvm/tree/sycl/sycl/test
