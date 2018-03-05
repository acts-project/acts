# A Common Tracking Software (ACTS) Project {#mainpage}

[![build status](https://gitlab.cern.ch/acts/acts-core/badges/master/build.svg)](https://gitlab.cern.ch/acts/acts-core/commits/master)
[![coverage report](https://gitlab.cern.ch/acts/acts-core/badges/master/coverage.svg)](https://gitlab.cern.ch/acts/acts-core/commits/master)

1. [Introduction](#intro)
2. [Mailing list](#mailing-list)
3. [Getting started](#getting-started)
    1. [Prerequisites](#prerequisites)
    2. [Installation](#installation)
    3. [CMake build system](#cmake)
    4. [Build ACTS on lxplus](#build-lxplus)
    5. [Build ACTS on your local machine](#build-local)
4. [Using ACTS in your own cmake project](#using-acts)
5. [Documentation](#documentation)
6. [License and authors](#license-authors)

# <a name="intro">Introduction</a>

This project is supposed to be an experiment-independent set of track reconstruction tools. The main philosophy is to provide high-level track reconstruction modules that can be used for any tracking detector. The description of the tracking detector's geometry is optimized for efficient navigation and quick extrapolation of tracks. Converters for several common geometry description languages exist. Having a highly performant, yet largely customizable implementation of track reconstruction algorithms was a primary objective for the design of this toolset. Additionally, the applicability to real-life HEP experiments played a major role in the development process. Apart from algorithmic code, this project also provides an event data model for the description of track parameters and measurements.

Key features of this project include:
* tracking geometry description which can be constructed from TGeo, DD4Hep, or gdml input,
* simple and efficient event data model,
* performant and highly flexible algorithms for track propagation and fitting,
* basic seed finding algorithms.

The git repository for the ACTS project can be found at <a href="https://gitlab.cern.ch/acts/acts-core.git">https://gitlab.cern.ch/acts/acts-core.git</a>.

# <a name="mailing-list">Mailing list</a>

In order to receive the latest updates, users of the ACTS project are encouraged to subscribe to [acts-users@cern.ch](https://e-groups.cern.ch/e-groups/Egroup.do?egroupName=acts-users). This list provides:
- regular updates on the software,
- access to the ACTS JIRA project for bug fixes/feature requests,
- a common place for asking any kind of questions.

If you find a bug, have a feature request, or want to contribute to ACTS, please have a look at the [contribution guide](CONTRIBUTING.md).

# <a name="getting-started">Getting started</a>

## <a name="prerequisites">Prerequisites</a>

Only few dependencies are required to build the Core library of ACTS. A list of
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

Compatible versions of all dependencies are provided by the [LCG91
Release](http://lcginfo.cern.ch/release/91). This release is also
used in the continous integration system to test the software. A setup script
is provided that can be used to setup the environment on lxplus machines at
CERN

    source CI/setup_lcg91.sh

## <a name="installation">Installation</a>

The ACTS repository is hosted on the GitLab instance at CERN. For the time being
you need to have a full CERN account in order to access the repository. We are
working on a solution for non-CERN users. In order to aquire the latest version
from the git repository you can follow the instructions below.

    git clone https://gitlab.cern.ch/acts/acts-core.git <ACTS_DIR>

## <a name="cmake">CMake build system</a>

CMake is used as build system for compiling and installing ACTS.
For a complete list of cmake options please refer to the [official documentation](https://cmake.org/cmake/help/v3.1/index.html) and this nice [list of general cmake options](https://cmake.org/Wiki/CMake_Useful_Variables).
Important options relevant for the ACTS project are given below. They are set by adding '-D\<OPTION\>=\<VALUE\>' to the cmake command.

|option|default|description|
|------|-------|-----------|
|ACTS_BUILD_LEGACY            | ON  | build Legacy package |
|ACTS_BUILD_DOC               | OFF | build documentation |
|ACTS_BUILD_EXAMPLES          | ON  | build examples |
|ACTS_BUILD_TESTS             | ON  | build unit tests |
|ACTS_BUILD_INTEGRATION_TESTS | OFF | build integration tests |
|ACTS_BUILD_DD4HEP_PLUGIN     | OFF | build DD4HEP plugins |
|ACTS_BUILD_MATERIAL_PLUGIN   | ON  | build Material plugins |
|ACTS_BUILD_TGEO_PLUGIN       | OFF | build TGeo plugins |
|CMAKE_INSTALL_PREFIX         |     | target installation directory |
|CMAKE_PREFIX_PATH            |     | search path for external packages |
|CMAKE_CXX_COMPILER           |     | set C++ compiler (e.g. g++ or clang++) |
|CMAKE_BUILD_TYPE             |     | build type (e.g. Debug, Release) affects compiler flags |
|DD4hep_DIR                   |     | path to the DD4hep installation |

## <a name="build-lxplus">Build ACTS on lxplus</a>

On lxplus the dependencies are provided by a LCG release. You can use the
following commands to build ACTS with all plugins using the same dependency
versions as in the continous integration system.

    source CI/setup_lcg89.sh
    mkdir build && cd build
    cmake -DCMAKE_INSTALL_PREFIX=<path you want> \
          -DACTS_BUILD_DD4HEP_PLUGIN=ON \
          -DACTS_BUILD_MATERIAL_PLUGIN=on \
          -DACTS_BUILD_TGEO_PLUGIN=ON ..
    make install

## <a name="build-local">Build ACTS on your local machine</a>

Building and running ACTS on your local machine is not offically supported.
However, if you have the necessary prerequisites installed it should be
possible to use it locally. ACTS developers regularly use different
recent Linux distributions and MacOS to build and develop ACTS.

# <a name="using-acts">Using ACTS in your own cmake project</a>

When using ACTS in your own cmake-based project, you need to include the following lines in your `CMakeLists.txt` file:

    find_package (ACTS COMPONENTS comp1 comp2 ...)

where `compX` are the required components from the ACTS project. See the `cmake` output for more information about which components are available.

# <a name="documentation">Documentation</a>

You can find a complete documentation of the ACTS functionality and the class reference guide at [http://acts.web.cern.ch/ACTS/latest/doc/index.html](http://acts.web.cern.ch/ACTS/latest/doc/index.html).

# <a name="license-authors">License and authors</a>

This project is published under the Mozilla Public License, v. 2.0. Details of
this license can be found in the [LICENSE](LICENSE) file or at
http://mozilla.org/MPL/2.0/.

Contributors to the ACTS project are listed in [AUTHORS](AUTHORS).

The ACTS project is based on the ATLAS tracking software. A list of contributors
to the ATLAS tracking repository can be found <a href="http://acts.web.cern.ch/ACTS/ATLAS_authors.html">here</a>.

The ACTS project contains a copy of [gcovr](http://gcovr.com) licensed under
the 3-Clause BSD license.
