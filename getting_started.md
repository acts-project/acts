# <a name="getting-started">Getting started</a>

1. [Getting started](#getting-started)
    1. [Prerequisites](#prerequisites)
    2. [Installation](#installation)
    3. [CMake build system](#cmake)
    4. [Build Acts on lxplus](#build-lxplus)
    5. [Build Acts on your local machine](#build-local)
2. [Using Acts in your own cmake project](#using-acts)
3. [Documentation](#documentation)

## <a name="prerequisites">Prerequisites</a>

Only few dependencies are required to build the Core library of Acts. A list of
prerequisites required is given below with the version numbers indicating which
versions were tested. Older versions may or may not work, feedback is very
welcome.

The following dependencies are required:

*   A C++14 compatible compiler, e.g. [gcc](https://gcc.gnu.org) (>= 6.2) or [clang](http://clang.llvm.org) (>= 4.0)
*   [cmake](https://cmake.org) (>= 3.7)
*   [boost](http://boost.org) (>= 1.62, with `program_options` and `unit_test_framework`)
*   [Eigen](http://eigen.tuxfamily.org) (>= 3.2.9)

The following dependencies are optional and are only needed for some of the
components:

*   [DD4Hep](http://dd4hep.cern.ch) (>= 1.02) for the DD4Hep plugin
*   [doxygen](http://doxygen.org) (>= 1.8.11) for the documentation
*   [graphviz](http://www.graphviz.org) (>= 2.28.00) for the documentation
*   [ROOT](https://root.cern.ch) (>= 6.10.00) for the TGeo plugin

Compatible versions of all dependencies are provided by **LCG releases**.
The current recommended release for building Acts is the
[LCG94 Release](http://lcginfo.cern.ch/release/94). This release is also used in the
continous integration (CI) system to test the software. Setup scripts are provided
in the repository that can be used to setup this release, and a few others, on
lxplus machines at CERN (see [below](#installation)).


## <a name="installation">Installation and quick start</a>

The Acts repository is hosted on the GitLab instance at CERN. In order to aquire the latest
version from the git repository you can simply clone:

```bash
git clone https://gitlab.cern.ch/acts/acts-core.git <ACTS_DIR>
```

You can then `cd <ACTS_DIR>` continue building Acts:

```bash
source CI/setup_lcg94.sh
mkdir build && cd build
cmake ..
cmake --build . -- install
```


## <a name="cmake">CMake build system</a>

CMake is used as build system for compiling and installing Acts.  For a
complete list of CMake options please refer to the [official documentation](https://cmake.org/cmake/help/v3.1/index.html) 
and this nice [list of general cmake options](https://cmake.org/Wiki/CMake_Useful_Variables).
Important options relevant for the Acts project are given below. They are set
by adding `-D<OPTION>=<VALUE>` to the `cmake` command.

| option                           | default | description                                             |
|----------------------------------|---------|---------------------------------------------------------|
| ACTS_BUILD_LEGACY                | ON      | build Legacy package                                    |
| ACTS_BUILD_DOC                   | OFF     | build documentation                                     |
| ACTS_BUILD_EXAMPLES              | OFF     | build examples                                          |
| ACTS_BUILD_TESTS                 | ON      | build unit tests                                        |
| ACTS_BUILD_INTEGRATION_TESTS     | OFF     | build integration tests                                 |
| ACTS_BUILD_DIGITIZATION_PLUGIN   | OFF     | build geometric digitization plugin                     |
| ACTS_BUILD_DD4HEP_PLUGIN         | OFF     | build DD4HEP plugin for DD4hep geometry                 |
| ACTS_BUILD_TGEO_PLUGIN           | OFF     | build TGeo plugin for ROOT geometry                     |
| ACTS_BUILD_JSON_PLUGIN           | OFF     | build Json plugin for Json geometry input/output        |
| ACTS_BUILD_MATERIAL_PLUGIN       | OFF     | build Material plugin                                   |
| ACTS_BUILD_IDENTIFICATION_PLUGIN | OFF     | build Identification plugin                             |
| CMAKE_INSTALL_PREFIX             |         | target installation directory                           |
| CMAKE_PREFIX_PATH                |         | search path for external packages                       |
| CMAKE_CXX_COMPILER               |         | set C++ compiler (e.g. g++ or clang++)                  |
| CMAKE_BUILD_TYPE                 |         | build type (e.g. Debug, Release) affects compiler flags |
| DD4hep_DIR                       |         | path to the DD4hep installation                         |

## <a name="building-acts">Building Acts</a>

### <a name="build-lxplus">Building Acts on lxplus</a>

The first step to build Acts is to acquire supported versions of the
dependencies.  Which dependencies are required depends on the plugins you
enable, as mentioned above.

If you are in an *lxplus-like* environment (i.e. `SLC6` or `CC7`, with
`cvmfs`), you can use the setup scripts located in `<ACTS_DIR>/CI` to get the
dependencies:


```bash
source <ACTS_DIR>/CI/setup_lcgXYZ.sh
```

where `XYZ` is the version number of the LCG release. There are multiple setup
scripts for different LCG releases, which corresponds to the releases the CI
tests against. The releases which can be set up using these scripts are therefore
sure to be compatible. **You can build Acts with any of these releases**. 
Additionally, there is a script called `setup_clang.sh` which will make the `clang` compiler available on top of one
of the LCG releases. This configuration is also tested by the CI.

Using one of the scripts, you can use the following commands to build Acts with
all plugins using the same dependency versions as in the continous integration
system.

```bash
source CI/setup_lcg94.sh # example, you can use any of the provided scripts.
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=<path you want> \
      -DACTS_BUILD_DD4HEP_PLUGIN=ON \
      -DACTS_BUILD_MATERIAL_PLUGIN=ON \
      -DACTS_BUILD_TGEO_PLUGIN=ON ..
cmake --build . -- install
```

In this example the DD4hep, Material and TGeo plugins. The install prefix is
set to `<path you want>`.

### <a name="build-local">Building Acts on your local machine</a>

Building and running Acts on your local machine is not offically supported.
However, if you have the necessary prerequisites installed it should be
possible to use it locally. Acts developers regularly use different
recent Linux distributions and macOS to build and develop Acts.

# <a name="using-acts">Using Acts in your own cmake project</a>

When using Acts in your own cmake-based project, you need to include the
following lines in your `CMakeLists.txt` file:

```bash
find_package (Acts COMPONENTS comp1 comp2 ...)
```

where `compX` are the required components from the Acts project. See the
`cmake` output for more information about which components are available.

# <a name="documentation">Documentation</a>

You can find a complete documentation of the Acts functionality and the class reference guide at [http://acts.web.cern.ch/ACTS/latest/doc/index.html](http://acts.web.cern.ch/ACTS/latest/doc/index.html).

