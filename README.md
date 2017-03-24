# A Common Tracking Software (ACTS) Project {#mainpage}

1. [Introduction](#intro)
2. [Mailing list](#mailing-list)
3. [Getting started](#getting-started)
    1. [Prerequisites](#prerequisites)
    2. [Installation](#installation)
    3. [cmake options](#cmake-options)
    4. [Using docker](#using-docker)
    5. [Building on lxplus](#lxplus-build)
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

+ [clang](http://clang.llvm.org/) (>= 3.8.0) or [gcc](https://gcc.gnu.org/) (>= 6.2)
+ [cmake](https://cmake.org/) (>= 3.5)
+ [boost](http://boost.org/) (>= 1.62, with <tt>program_options</tt>)
+ [Eigen](http://eigen.tuxfamily.org/) (>= 3.2.10)

The following dependencies are optional and are only needed for some of the plugins

+ [doxygen](http://doxygen.org) (>= 1.8.11) for the documentation
+ [graphviz](http://www.graphviz.org/) (>= 2.26.00) for the documentation
+ [ROOT](https://root.cern.ch/) (>= 6.08.00) for TGeo plugin & for DD4hep plugin
+ [DD4Hep](https://github.com/AIDASoft/DD4hep) for DD4Hep plugin

### Using recent CMake releases from CVMFS

Recent releases of CMake can be found on CVMFS at `/cvmfs/sft.cern.ch/lcg/contrib/CMake/`. These builds are self-contained and intended to work even on non-Redhat linux distributions. In particular, they have been successfully used on Ubuntu 14.04.

To use the CMake release x.y.z from this source, do...

    export PATH=/cvmfs/sft.cern.ch/lcg/contrib/CMake/x.y.z/Linux-x86_64/bin:${PATH}

## <a name="installation">Installation</a>

The ACTS repository is hosted on the GitLab instance at CERN. For the time being you need to have a full CERN account in order to access the repository. We are working on a solution for non-CERN users. Cmake is used as build system for compiling and installing ACTS libraries. For a complete description of the available cmake options please see below.

In order to install the latest version, you can follow the instructions below where \<DIR\> refers to some directory which needs to be set depending on your local system configuration.

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
|BOOST_ROOT             | empty                 | path to the ROOT installation                           |
|EIGEN_INCLUDE_DIR      | empty                 | path to the Eigen installation                          |
|BUILD_DD4HEP_PLUGIN    | OFF                   | build DD4Hep plugin (requires TGeoPlugin, ROOT, DD4hep) |
|BUILD_DOC              | OFF                   | build doxygen documentation (requires doxygen)          |
|BUILD_TESTS            | ON                    | build unit tests                                        |
|BUILD_TGEO_PLUGIN      | OFF                   | build TGeo plugin (requires ROOT)                       |
|BUILD_MATERIAL_PLUGIN  | OFF                   | build material mapping plugin                           |
|CMAKE_INSTALL_PREFIX   | empty                 | target installation directory                           |
|CMAKE_PREFIX_PATH      | empty                 | search path for external packages                       |
|CMAKE_CXX_COMPILER     | empty                 | set C++ compiler (e.g. g++ or clang++)                  |
|CMAKE_BUILD_TYPE       | None                  | build type (e.g. Debug, Release) affects compiler flags |

## <a name="using-docker">Build ACTS using docker</a>

The ACTS team provides you with a [docker](https://en.wikipedia.org/wiki/Docker_(software)) image with all required software already pre-installed. This is the very same image used in our continuous integration system. Hence, it is very well tested and should be the easiest way to get you started. In order to use it, you need to have docker installed. On Ubuntu systems one could achieve this by running

> sudo apt-get install docker.io

While the docker image provides you with the environment for building ACTS, it does not contain the source code itself. The reasoning behind this is that you can develop ACTS on your host machine using your preferred development tools/editors/GUIs and use the docker container only for compiling/testing. Therefore, you need to clone the ACTS repository first

> git clone https://gitlab.cern.ch/acts/a-common-tracking-sw.git acts

As a second step you need to pull the ACTS docker image

> docker pull gitlab-registry.cern.ch/acts/a-common-tracking-sw

Before starting the docker container, you can create a shorter tag for this image to avoid a lot of typing

> docker tag gitlab-registry.cern.ch/acts/a-common-tracking-sw acts

Now spin up the docker container with the mysterious command

> docker run -d -t -i -v acts:/acts -e LOCAL_USER_ID=\`id -u\` -e LOCAL_GROUP_ID=\`id -g\` --name acts acts

Here is what it means:

- -d runs the container in the background (detached state)
- -t gives you acces to a shell (bash)
- -i stands for interactive and allows you to attach
- -v maps the directory `acts` from your host machine to `/acts` inside the container
- -e sets some environment variables which are used to map the current user to the user inside the container
- --name gives a name to the container
- the last argument is a reference to the docker image used for creating this container

You can attach to the container using

> docker attach acts

You can then go ahead like you would on your host machine and start building ACTS using `cmake ... && make`. Remember that the ACTS source code is located under `/acts` inside the container. There is also a simple python wrapper script called `acts-build` in case you do not remember the longish cmake command. Running `acts-build --help` gives you a (short) list of available options.  

For instance you could compile and install ACTS using

> acts-build /acts/ /workdir/build --make-options "install" --cmake-options " -DCMAKE_INSTALL_PREFIX=/workdir/install"<br />
> cd /workdir/build && make test<br />
> cd /workdir/install/bin<br />
> source setup.sh<br />
> ./Examples/ACTSGenericDetector

You can detach from the container again pressing the key sequence `CTRL+P CTRL+Q`.  
If you just want to test the compilation non-interactively, you could also execute (from the host machine)

> docker exec acts acts-start acts-build /acts /workdir/build

This command could, for instance, be used as custom build command in IDEs.

## <a name="lxplus-build">Build ACTS on lxplus</a>

On lxplus many of the dependencies are provided by the LCG releases. A possible setup is the following:

> export actsbase=${PWD}<br />
> export lcgversion=LCG_87<br />
> export platform=x86_64-slc6-gcc62-opt<br />
> export lcgdir=/cvmfs/sft.cern.ch/lcg/views/${lcgversion}/${platform}<br />
> source /cvmfs/sft.cern.ch/lcg/views/${lcgversion}/${platform}/setup.sh<br />
> cd /cvmfs/sft.cern.ch/lcg/releases/${lcgversion}/DD4hep/00-17/${platform} && source bin/thisdd4hep.sh && cd ${actsbase}<br />
> <br />
> mkdir build && cd build<br />
> cmake .. -DCMAKE_INSTALL_PREFIX=<path you want> \<br />
>	 -DEIGEN_INCLUDE_DIR=${lcgdir}/include/eigen3/ \<br />
>	-DBOOST_INCLUDEDIR=${lcgdir}/include/boost-1_62/ \<br />
>	-DBUILD_TGEO_PLUGIN=ON \<br />
>	-DBUILD_DD4HEP_PLUGIN=ON<br />
> make install

# <a name="using-acts">Using ACTS in your own cmake project</a>

When using ACTS in your own cmake-based project, you need to include the following lines in your `CMakeLists.txt` file:

> find_package (ACTS COMPONENTS comp1 comp2 ...)

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
