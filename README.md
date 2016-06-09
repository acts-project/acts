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

+ [clang](http://clang.llvm.org/) (>= 3.8.0) or [gcc](https://gcc.gnu.org/) (>= 5.3.1)
+ [cmake](https://cmake.org/) (>= 2.8)
+ [boost](http://www.boost.org/) (>= 1.61)
+ [Eigen](http://eigen.tuxfamily.org/) (>= 3.2.8)

The following dependencies are optional and only needed for some of the plugins

+ [ROOT](https://root.cern.ch/) (>= 6.06.04)
+ DD4Hep

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

For a complete list of cmake options please refer to the [official documentation](https://cmake.org/cmake/help/v3.1/index.html).
A full list of ACTS specific cmake options can be obtained by running the following command
> cmake \<ACTS_DIR\> -DPRINT_OPTIONS=ON
Important options relevant for the ACTS project are given below. They can be set
by adding '-D<OPTION>=<VALUE>' to the cmake command.

|option|default|description|
|------|-------|-----------|
|BOOST_ROOT| empty | |


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

