[![Build Status](http://acts-ci:8080/job/ACTS-CI/badge/icon)](http://acts-ci:8080/job/ACTS-CI)

# A Common Tracking Software (ACTS) Project

This library is based on the track reconstruction software developed by the 
[ATLAS Collaboration](http://cern.ch/atlas).

The main philosophy is to provide high-level track reconstruction modules that 
can be specified for detector technologies by simple extension.

* event data model (EDM)
* geometry 
* track reconstruction tools

The library is attempted to build against Gaudi and Gaudi-Athena, while
additional external dependencies are kept at a minimum.


Please visit our [Wiki](https://gitlab.cern.ch/acts/a-common-tracking-sw/wikis/home) for more information.

# Getting Started

## Prerequisites

Only few dependencies are required to build the Core library of ACTS. The given
version numbers were tested to be working. Older versions may or may not work.

The following dependencies are required:

+ [clang](http://clang.llvm.org/) (>= 3.8.0) or [gcc](https://gcc.gnu.org/) (>= 4.9.3)
+ [cmake](https://cmake.org/) (>= 2.8)
+ [boost](http://boost.org/) (>= 1.59)
+ [Eigen](http://eigen.tuxfamily.org/) (>= 3.2.8)

The following dependencies are optional and only needed for some of the plugins

+ [doxygen](http://doxygen.org) (>= 1.6.1) for the documentation
+ [graphviz](http://www.graphviz.org/) (>= 2.26.00) for the documentation
+ [ROOT](https://root.cern.ch/) (>= 6.06.04) for TGeo plugin
+ DD4Hep for DD4Hep plugin

## Installation

The ACTS repository is hosted on the GitLab instance at CERN. For the time being
you need to have a full CERN account in order to access this. We are working on
a solution for non-CERN users. We use cmake as build system for compiling and
installing ACTS libraries. For a complete description of the available cmake
options please see below.

In order to install the latest version, you can follow the instructions below
where \<DIR\> refers to some directory which needs to be set.

> git clone https://gitlab.cern.ch/acts/a-common-tracking-sw.git \<ACTS_DIR\><br />
> mkdir \<BUILD_DIR\><br />
> cd \<BUILD_DIR\><br />
> cmake .. -DEIGEN_INCLUDE_DIR=\<EIGEN_INSTALLATION\> -DBOOST_ROOT=\<BOOST_INSTALLATION\><br />
> make<br />
> make install<br />

## cmake options

For a complete list of cmake options please refer to the [official documentation](https://cmake.org/cmake/help/v3.1/index.html)
and this nice [list of general cmake options](https://cmake.org/Wiki/CMake_Useful_Variables).
A full list of ACTS specific cmake options can be obtained by running the following command
> cmake \<ACTS_DIR\> -DPRINT_OPTIONS=ON
Important options relevant for the ACTS project are given below. They can be set
by adding '-D<OPTION>=<VALUE>' to the cmake command.

|option|default|description|
|------|-------|-----------|
|BOOST_ROOT           | empty                 | path to the ROOT installation          |
|EIGEN_INCLUDE_DIR    | empty                 | path to the Eigen installation         |
|BUILD_DD4HEP_PLUGIN  | OFF                   | build DD4Hep plugin (requires DD4hep)  |
|BUILD_DOC            | OFF                   | build doxygen documentation            |
|BUILD_TESTS          | ON                    | build unit tests                       |
|BUILD_TGEO_PLUGIN    | OFF                   | build TGeo plugin (requires ROOT)      |
|CMAKE_INSTALL_PREFIX | /opt/ACTS/\<version\> | target installation directory          |
|CMAKE_PREFIX_PATH    | empty                 | search path for external packages      |    
|CMAKE_CXX_COMPILER   | empty                 | set C++ compiler (e.g. g++ or clang++) |    
|CMAKE_CXX_FLAGS      | -std=c++14            | set C++ compiler flags (e.g. "-O2 -g") |  
|CMAKE_BUILD_TYPE     | None                  | build type (e.g. Debug, Release)       |

## Example build on lxplus at CERN

If you are logged in to lxplus at CERN, you could run the following commands
to install ACTS with all its components and run an example. Please note the 
ong cmake command which spans several lines.

> git clone ssh://git@gitlab.cern.ch:7999/acts/a-common-tracking-sw.git acts<br />
> mkdir acts/build<br />
> cd acts/build<br />
> source /afs/cern.ch/sw/lcg/contrib/gcc/4.9.3/x86_64-slc6/setup.sh<br />
> cmake .. \\ <br />
>   -DEIGEN_INCLUDE_DIR=/afs/cern.ch/sw/lcg/releases/eigen/3.2.7-292e1/x86_64-slc6-gcc49-opt/include/eigen3/ \\ <br />
>   -DBOOST_ROOT=/afs/cern.ch/sw/lcg/releases/LCG_83/Boost/1.59.0_python2.7/x86_64-slc6-gcc49-opt/include/boost-1_59/ \\ <br />
>   -DBUILD_DOC=ON \\ <br />
>   -DBUILD_TGEO_PLUGIN=ON \\ <br />
>   -DBUILD_DD4HEP_PLUGIN=ON \\ <br />
>   -DCMAKE_PREFIX_PATH="/afs/cern.ch/exp/fcc/sw/0.7/DD4hep/20161003/x86_64-slc6-gcc49-opt/;/afs/cern.ch/sw/lcg/releases/LCG_83/ROOT/6.06.00/x86_64-slc6-gcc49-opt/cmake" \\ <br />
>   -DCMAKE_INSTALL_PREFIX=\`pwd\`/installed <br />
> make<br />
> make doc<br />
> make install<br />
> cd installed <br />
> export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:\`pwd\`/lib64 <br />
> ./bin/ACTSGenericDetector

# License and Authors

This project is published under the Mozilla Public License, v. 2.0. Details of
this license can be found in the LICENSE file or at http://mozilla.org/MPL/2.0/.

The authors of this project are (in alphabetical order):
- Noemi Calace
- Christian Gumpert
- Julia Hrdinka
- Andreas Salzburger

Some of the code was originally written for the ATLAS software. A list of
contributors for this portion of the code will be added shortly.

