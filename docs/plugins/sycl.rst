SYCL plugin
===========

The SYCL plugin allows building interfaces that implement algorithms in SYCL.

Build
-----

During the CMake configuration the flag ``CMAKE_BUILD_PLUGIN_SYCL=on`` needs to be set. CMake will check whether our compiler is SYCL compatible and will look for available SYCL target architectures. By default, SPIR 64-bits targets are looked for, but this can be configured with the ``SYCL_POSSIBLE_TARGETS`` variable. For example, to target CUDA backends, we should set this variable to ``nvptx64-nvidia-cuda-sycldevice``.

Also, by setting the variable ``SYCL_SEARCH_SUPPORT_LIBRARIES`` it is possible to specify support libraries, or object files that are needed for execution on a specific target. Missing object files cause a runtime error.

Device selection
----------------

Fisrtly, we can list available platforms and devices with **clinfo**:

.. code-block:: bash

    [bash][atspot01]:build > clinfo
        Number of platforms                             3
        Platform Name                                   Intel(R) OpenCL
        Platform Vendor                                 Intel(R) Corporation
        Platform Version                                OpenCL 2.1 LINUX
        ...

CMake should recognize these platforms during configuration (if not told otherwise).

There are more options we can choose from to configure device selection (meaning constructing a ``cl::sycl::device_selector`` object).

We can alter the behavior of the SYCL ``default_selector`` by setting the environment variable ``SYCL_BE`` (BE stands for backend). For example by setting it to ``PI_CUDA`` it forces the the usage of the CUDA backend (if available). Similarly ``PI_OPENCL`` forces the usage of the OpenCL backend. For more information see this Getting Started `guide`_.:

.. code-block:: bash

    SYCL_BE=PI_CUDA <binary> <arguments>

We can also use a custom device selector. As an example, see :class:`Acts::Sycl::DeviceSelector`, which selects a non OpenCL CUDA backend. This is the device selector used by the :class:`Acts::Sycl::QueueWrapper` object that manages the construction, ownership and destruction of a ``cl::sycl::queue`` object. This solution allows us to construct a ``queue`` object in translation units that are not linked against SYCL.

.. _guide: https://intel.github.io/llvm-docs/GetStartedGuide.html#run-simple-dpc-application

Algorithms
----------
Currently, the following algorithms are implemented in SYCL:
 #. Seeding

SYCL Seeding
------------
**Usage**

To create a Seedfinder object that implements the seeding algorithm in SYCL, we need to instantiate an object from :class:`Acts::Sycl::Seedfinder`. Similarly to :class:`Acts::Seedfinder`, we need to provide a configuration object of :class:`Acts::SeedfinderConfig<external_spacepoint_t>` and an object that implements *experiment specific* cuts. As these cuts are performed in SYCL kernels (on the device side), instead of a :class:`IExperimentCuts` instance, we need to construct a :class:`DeviceExperimentCuts` one.

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

Optionally we can also give a :class:`Acts::Sycl::QueueWrapper` object to the constructor of :class:`Acts::Sycl::Seedfinder`, which is a wrapper object around a ``cl::sycl::queue`` type. This allows us to construct our own ``queue`` instance and to reuse it.

**Performance**

It depends on our architecture, the size of the event we are reconstructing and the effectiveness of our experiment specific cuts how well the algorithm performs, and whether we can benefit at all from using the SYCL plugin. It is advised to compare perfomance and precision first with the CPU version of the choosen algorithm,. This should be possible with the tests provided.

Resources
---------

For more information about SYCL see the `specification`_ (date: 2020. September 7.).
There is a `documentation`_ for Intel implementation and DPC++ extensions (see `examples`_).

.. _documentation: https://software.intel.com/content/www/us/en/develop/download/intel-oneapi-programming-guide.html
.. _specification: https://www.khronos.org/registry/SYCL/specs/sycl-1.2.1.pdf
.. _examples: https://github.com/intel/llvm/tree/sycl/sycl/test
