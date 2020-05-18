# Getting started

## Prerequisites

The following dependencies are required to build the Acts Core library:

*   A C++17 compatible compiler (recent versions of either gcc and clang should work)
*   [CMake](https://cmake.org) (>= 3.11)
*   [Boost](http://boost.org) (>= 1.69, with `filesystem`, `program_options`, and `unit_test_framework`)
*   [Eigen](http://eigen.tuxfamily.org) (>= 3.2.9)

The following dependencies are optional and are needed to build additional
components.

*   [DD4Hep](http://dd4hep.cern.ch) (>= 1.10) for the DD4Hep plugin and some examples
*   [Doxygen](http://doxygen.org) (>= 1.8.11) for the documentation
*   [Geant4](http://geant4.org/) for some examples
*   [HepMC](https://gitlab.cern.ch/hepmc/HepMC3) (>= 3.1) for some examples
*   [Pythia8](http://home.thep.lu.se/~torbjorn/Pythia.html) for some examples
*   [Intel Threading Building Blocks](https://01.org/tbb) for the examples
*   [ROOT](https://root.cern.ch) (>= 6.10) for the TGeo plugin and the examples
*   [Sphinx](https://www.sphinx-doc.org) (>= 2.0) for the documentation

Compatible versions of all dependencies are provided by **LCG releases**.
The current recommended release for building Acts is
[LCG 95](http://lcginfo.cern.ch/release/96). This release is also used in the
continous integration (CI) system to test the software. Setup scripts are provided
in the repository that can be used to setup this release, and a few others, on
lxplus machines at CERN (see [below](#installation)).

## Quick start

The Acts repository is hosted Github. In order to aquire the latest
version from the git repository you can simply clone:

```console
git clone https://github.com/acts-project/acts.git
```

You can then `cd acts` to continue building Acts:

```console
source CI/setup_lcg94.sh
mkdir build && cd build
cmake ..
cmake --build . -- install
```

## CMake build system

CMake is used as build system for compiling and installing Acts. For a
complete list of CMake options please refer to the [official documentation](https://cmake.org/cmake/help/v3.1/index.html)
and this nice [list of general cmake options](https://cmake.org/Wiki/CMake_Useful_Variables).
Important options relevant for the Acts project are given below. They are set
by adding `-D<OPTION>=<VALUE>` to the `cmake` command.

| option                           | default | description                                             |
|----------------------------------|---------|---------------------------------------------------------|
| ACTS_BUILD_DD4HEP_PLUGIN         | OFF     | Build DD4HEP geometry plugin                            |
| ACTS_BUILD_DIGITIZATION_PLUGIN   | OFF     | Build Digitization plugin                               |
| ACTS_BUILD_JSON_PLUGIN           | OFF     | Build Json plugin for geometry input/output             |
| ACTS_BUILD_IDENTIFICATION_PLUGIN | OFF     | Build Identification plugin                             |
| ACTS_BUILD_TGEO_PLUGIN           | OFF     | Build TGeo geometry plugin                              |
| ACTS_BUILD_FATRAS                | OFF     | Build FAst TRAcking Simulation package                  |
| ACTS_BUILD_LEGACY                | OFF     | Build Legacy package                                    |
| ACTS_BUILD_BENCHMARKS            | OFF     | Build benchmarks                                        |
| ACTS_BUILD_EXAMPLES              | OFF     | Build examples                                          |
| ACTS_BUILD_UNITTESTS             | OFF     | Build unit tests                                        |
| ACTS_BUILD_INTEGRATIONTESTS      | OFF     | Build integration tests                                 |
| ACTS_BUILD_DOC                   | OFF     | Build documentation                                     |
| ACTS_USE_BUNDLED_NLOHMANN_JSON   | ON      | Use external or bundled Json library                    |
| CMAKE_INSTALL_PREFIX             |         | The installation directory                              |
| CMAKE_PREFIX_PATH                |         | Search path for external packages                       |
| CMAKE_CXX_COMPILER               |         | Set C++ compiler (e.g. g++ or clang++)                  |
| CMAKE_BUILD_TYPE                 |         | Build type (e.g. Debug, Release) affects compiler flags |
| DD4hep_DIR                       |         | Path to the DD4hep installation                         |

## Building Acts

### On lxplus

The first step to build Acts is to acquire supported versions of the
dependencies.  Which dependencies are required depends on the plugins you
enable, as mentioned above.

If you are in an *lxplus-like* environment (i.e. `SLC6` or `CC7`, with
`cvmfs`), you can use the setup scripts located in `<ACTS_DIR>/CI` to get the
dependencies:

```console
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

```console
source CI/setup_lcg94.sh # example, you can use any of the provided scripts.
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=<path you want> \
      -DACTS_BUILD_DD4HEP_PLUGIN=ON \
      -DACTS_BUILD_TGEO_PLUGIN=ON ..
cmake --build . -- install
```

In this example the DD4hep, Material and TGeo plugins. The install prefix is
set to `<path you want>`.

### On your local machine

Building and running Acts on your local machine is not offically supported.
However, if you have the necessary prerequisites installed it should be
possible to use it locally. Acts developers regularly use different
recent Linux distributions and macOS to build and develop Acts.

## Using Acts

When using Acts in your own CMake-based project, you need to include the
following lines in your `CMakeLists.txt` file:

```cmake
find_package (Acts COMPONENTS comp1 comp2 ...)
```

where `compX` are the required components from the Acts project. See the
`cmake` output for more information about which components are available.
