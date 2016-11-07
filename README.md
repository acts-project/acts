# A Common Tracking Software (ACTS) Project {#mainpage}

1. [Introduction](#intro)
2. [Mailing list](#mailing-list)
3. [Getting started](#getting-started)
    1. [Prerequisites](#prerequisites)
    2. [Installation](#installation)
    3. [cmake options](#cmake-options)
    4. [Example build on lxplus at CERN](#example-build)
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

The git repository for the ACTS project can be found at <a href="https://gitlab.cern.ch/acts/a-common-tracking-sw.git">https://gitlab.cern.ch/acts/a-common-tracking-sw.git</a>.

# <a name="mailing-list">Mailing list</a>

In order to receive the latest updates, users of the ACTS project are encouraged to subscribe to [acts-users@cern.ch](https://e-groups.cern.ch/e-groups/Egroup.do?egroupName=acts-users). This list provides:
- regular updates on the software,
- access to the ACTS JIRA project for bug fixes/feature requests,
- a common place for asking any kind of questions.

If you find a bug, have a feature request, or want to contribute to ACTS, please have a look at the [contribution guide](CONTRIBUTING.md).

# <a name="getting-started">Getting started</a>

## <a name="prerequisites">Prerequisites</a>

Only few dependencies are required to build the Core library of ACTS. A list of prerequisites required is given below with the version numbers indicating which versions were tested. Older versions may or may not work, feedback is very welcome.

The following dependencies are required:

+ [clang](http://clang.llvm.org/) (>= 3.8.0) or [gcc](https://gcc.gnu.org/) (>= 4.9.3)
+ [cmake](https://cmake.org/) (>= 2.8)
+ [boost](http://boost.org/) (>= 1.60)
+ [Eigen](http://eigen.tuxfamily.org/) (>= 3.2.8)

The following dependencies are optional and are only needed for some of the plugins

+ [doxygen](http://doxygen.org) (>= 1.8.11) for the documentation
+ [graphviz](http://www.graphviz.org/) (>= 2.26.00) for the documentation
+ [ROOT](https://root.cern.ch/) (>= 6.06.04) for TGeo plugin & for DD4hep plugin
+ [DD4Hep](https://github.com/AIDASoft/DD4hep) for DD4Hep plugin

## <a name="installation">Installation</a>

The ACTS repository is hosted on the GitLab instance at CERN. For the time being you need to have a full CERN account in order to access the repository. We are working on a solution for non-CERN users. Cmake is used as build system for compiling and installing ACTS libraries. For a complete description of the available cmake options please see below.

In order to install the latest version, you can follow the instructions below where \<DIR\> refers to some directory which needs to be set depending on your configuration.

> git clone https://gitlab.cern.ch/acts/a-common-tracking-sw.git \<ACTS_DIR\><br />
> mkdir \<BUILD_DIR\><br />
> cd \<BUILD_DIR\><br />
> cmake \<ACTS_DIR\> -DEIGEN_INCLUDE_DIR=\<EIGEN_INSTALLATION\> -DBOOST_ROOT=\<BOOST_INSTALLATION\> -DCMAKE_INSTALL_PREFIX=\<INSTALL_DIR\><br />
> make<br />
> make install<br />

## <a name="cmake-options">cmake options</a>

For a complete list of cmake options please refer to the [official documentation](https://cmake.org/cmake/help/v3.1/index.html) and this nice [list of general cmake options](https://cmake.org/Wiki/CMake_Useful_Variables). A full list of ACTS specific cmake options can be obtained by running the following command

> cmake \<ACTS_DIR\> -DPRINT_OPTIONS=ON

Important options relevant for the ACTS project are given below. They can be set by adding '-D\<OPTION\>=\<VALUE\>' to the cmake command.

|option|default|description|
|------|-------|-----------|
|BOOST_ROOT             | empty                 | path to the ROOT installation                  |
|EIGEN_INCLUDE_DIR      | empty                 | path to the Eigen installation                 |
|BUILD_DD4HEP_PLUGIN    | OFF                   | build DD4Hep plugin (requires DD4hep)          |
|BUILD_DOC              | OFF                   | build doxygen documentation (requires doxygen) |
|BUILD_TESTS            | ON                    | build unit tests                               |
|BUILD_TGEO_PLUGIN      | OFF                   | build TGeo plugin (requires ROOT)              |
|BUILD_MATERIAL_PLUGIN  | OFF                   | build material mapping plugin                  |
|CMAKE_INSTALL_PREFIX   | empty                 | target installation directory                  |
|CMAKE_PREFIX_PATH      | empty                 | search path for external packages              |     
|CMAKE_CXX_COMPILER     | empty                 | set C++ compiler (e.g. g++ or clang++)         |    
|CMAKE_BUILD_TYPE       | None                  | build type (e.g. Debug, Release) affects compiler flags |

## <a name="example-build">Example build on lxplus at CERN</a>

If you are logged in to lxplus at CERN, you can run the following commands to install the core components of ACTS and run an example. Please note the long cmake command which spans several lines.

> git clone ssh://git@gitlab.cern.ch:7999/acts/a-common-tracking-sw.git acts<br />
> mkdir acts/build<br />
> cd acts/build<br />
> source /afs/cern.ch/sw/lcg/contrib/gcc/4.9.3/x86_64-slc6/setup.sh<br />
> cmake .. \\ <br />
>   -DEIGEN_INCLUDE_DIR=/afs/cern.ch/sw/lcg/releases/eigen/3.2.7-292e1/x86_64-slc6-gcc49-opt/include/eigen3/ \\ <br />
>   -DBOOST_ROOT=/afs/cern.ch/sw/lcg/releases/LCG_83/Boost/1.59.0_python2.7/x86_64-slc6-gcc49-opt/include/boost-1_59/ \\ <br />
>   -DBUILD_DOC=ON \\ <br />
>   -DCMAKE_PREFIX_PATH="/afs/cern.ch/exp/fcc/sw/0.7/DD4hep/20161003/x86_64-slc6-gcc49-opt/;/afs/cern.ch/sw/lcg/releases/LCG_83/ROOT/6.06.00/x86_64-slc6-gcc49-opt/cmake" \\ <br />
>   -DCMAKE_INSTALL_PREFIX=\`pwd\`/installed <br />
> make<br />
> make doc<br />
> make install<br />
> cd installed <br />
> source bin/setup.sh <br />
> ./bin/ACTSGenericDetector

# <a name="using-acts">Using ACTS in your own cmake project</a>

When using ACTS in your own cmake-based project, you need to include the following lines in your `CMakeLists.txt` file:

> find_package (ACTS COMPONENTS comp1 comp2 ...)

where `compX` are the required components from the ACTS project.

# <a name="documentation">Documentation</a>

You can find a complete documentation of the ACTS functionality and the class reference guide at [http://acts.web.cern.ch/ACTS/latest/doc/index.html](http://acts.web.cern.ch/ACTS/latest/doc/index.html).

# <a name="license-authors">License and authors</a>

This project is published under the Mozilla Public License, v. 2.0. Details of
this license can be found in the [LICENSE](LICENSE) file or at
http://mozilla.org/MPL/2.0/.

Contributors to the ACTS project are listed in [AUTHORS](AUTHORS).

The ACTS project is based on the ATLAS tracking software. A list of contributors
to the ATLAS tracking repository can be found <a href="http://acts.web.cern.ch/ACTS/ATLAS_authors.html">here</a>.
