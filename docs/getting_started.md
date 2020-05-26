# Getting started

## Quick start

Acts is developed in C++ and is built using [CMake](https://cmake.org). Building
the core library requires a C++17 compatible compiler,
[Boost](http://boost.org), and [Eigen](http://eigen.tuxfamily.org). The
following commands will clone the repository, configure, and build the core
library:

```console
$ git clone https://github.com/acts-project/acts <source-dir>
$ cmake -B <build-dir> -S <source-dir>
$ cmake --build <build-dir>
```

For a full list of dependencies, including specific versions, see the
[Prerequisites](#prerequisites) section below. Build options to activate
additional components are described in the [Build options](#build-options)
section.

## Prerequisites

The following dependencies are required to build the Acts core library:

-   A C++17 compatible compiler (recent versions of either gcc and clang should work)
-   [CMake](https://cmake.org) >= 3.11
-   [Boost](http://boost.org) >= 1.69 with `filesystem`, `program_options`, and `unit_test_framework`
-   [Eigen](http://eigen.tuxfamily.org) >= 3.2.9

The following dependencies are optional and are needed to build additional
components.

-   [DD4Hep](http://dd4hep.cern.ch) >= 1.10 for the DD4Hep plugin and some examples
-   [Doxygen](http://doxygen.org) >= 1.8.11 for the documentation
-   [Geant4](http://geant4.org/) for some examples
-   [HepMC](https://gitlab.cern.ch/hepmc/HepMC3) >= 3.1 for some examples
-   [HepPDT](http://lcginfo.cern.ch/pkg/HepPDT) >= 2.06 for some examples
-   [Pythia8](http://home.thep.lu.se/~torbjorn/Pythia.html) for some examples
-   [Intel Threading Building Blocks](https://01.org/tbb) for the examples
-   [ROOT](https://root.cern.ch) >= 6.10 for the TGeo plugin and the examples
-   [Sphinx](https://www.sphinx-doc.org) >= 2.0 with [Breathe](https://breathe.readthedocs.io/en/latest/), [Exhale](https://exhale.readthedocs.io/en/latest/), and [recommonmark](https://recommonmark.readthedocs.io/en/latest/index.html) extensions for the documentation

All external dependencies must be provided prior to building Acts. Compatible
versions of all dependencies are provided e.g. by the [LCG
releases](http://lcginfo.cern.ch/). The minimum support release is [LCG
95apython3](http://lcginfo.cern.ch/release/95apython3). Other options are also
available and are discussed in the [Building Acts](#building-acts) section.

## Building Acts

Acts uses [CMake](https://cmake.org) to configure, build, and install the
software. After checking out the repository code into a `<source>` directory,
CMake is called first to configure the build into a separate `<build>`
directory. A typical setup is to create a `<source>/build` directory within the
sources, but this is just a convention; not a requirement. The following command
runs the configuration and searches for the dependencies. The `<build>`
directory is automatically created.

```console
$ cmake -B <build> -S <source>
```

The build can be configured via various options that are listed in detail in the
[Build options](#build-options) section. Options are set on the command line.
The previous command could be e.g. modified to

```console
$ cmake -B <build> -S <source> -DACTS_BUILD_UNITTESTS=on -DACTS_BUILD_FATRAS=on
```

After the configuration succeeded, the software is build. This is also done with cmake via the following command

```console
$ cmake --build <build>
```

This automatically calls the configure build tool, e.g. Make or Ninja. To build only a specific target, the target names has to be separated from the CMake options by `--`, i.e.

```console
$ cmake --build <build> -- ActsFatras # to build the Fatras library
```

The build commands are the same regardless of where you are building the
software. Depending on your build environment, there are different ways how to
make the dependencies available.

### On lxplus

On CERNs lxplus login machines, dependencies are provided via LCG releases
that are available through CVMFS. A setup script is provided to activate one of the compatible releases that can be used as follows:

```console
$ cd <source>
$ source CI/setup_lcg95.sh # or setup_lcg96.sh
```

This activates a compatible release variant regardless of whether you are on
*SLC 6* or on *CentOS 7*.

After sourcing the setup script, you can build Acts as described above. The
following commands will build Acts in the `<source>/build` directory with the
Fatras component.

```console
$ cd <source>
$ source CI/setup_lcg96.sh
$ cmake -B build -S . -DACTS_BUILD_FATRAS=on
$ cmake --build build
```

### In a container

A set of container images is available through the [Acts container
registry][acts_containers]. The following containers are used as part of the
continous integration setup and come with all dependencies pre-installed.

-   `ubuntu2004`: based on Ubuntu 20.04 with manual installation of HEP-specific
    software
-   `centos7-lcg95apython3`: based on CentOS 7 with HEP-specific software from
    LCG release 95apython
-   `centos7-lcg96`: based on CentOS 7 with HEP-specific software from LCG
    release 96
-   `slc6-lcg95apython3`: based on SLC 6 with HEP-specific software from LCG
    release 95apython

To use these locally, you first need to pull the relevant images from the
registry. The current version is always tagged as `master.

```console
$ docker pull gitlab-registry.cern.ch/acts/machines/ubuntu2004:master
```

This should print the image id as part of the output. You can also find out the
image id by running `docker images` to list all your locally available container
images. If you encounter authorization issues, you might need to login to the
registry first via

```console
$ docker login gitlab-registry.cern.ch
```

Now, you need to start a shell within the container to run the build. Assuming
that `<source>` is the path to your source checkout on your host machine, the
following command will make the source directory available as `/acts` in the
container and start an interactive `bash` shell

```console
$ docker run --volume=<source>:/acts:ro --interactive --tty <image> /bin/bash
```

where `<image>` is the image id that was previously mentioned. If you are using the Ubuntu-based image you are already good to go. For the images based on LCG releases, you can now activate the LCG release in the container shell by sourcing a setup script:

```console
container$ source /opt/lcg_view/setup.sh
```

Building Acts follows the instructions above with `/acts` as the source directory, e.g.

```console
container$ cmake -B build -S /acts -DACTS_BUILD_FATRAS=on
container$ cmake --build build
```

[acts_containers]: https://gitlab.cern.ch/acts/machines/container_registry


### On your local machine

Building and running Acts on your local machine is not offically supported.
However, if you have the necessary prerequisites installed it is possible to use
it locally. Acts developers regularly use different recent Linux distributions
and macOS to build and develop Acts.

## Building the documentation

The documentation uses [Doxygen][doxygen] to extract the source code
documentation and [Sphinx][sphinx] with the [Breathe][breathe] and
[Exhale][exhale] extensions to generate the documentation website. To build the
documentation locally, you need to have [Doxygen][doxygen] installed from your
package manager. [Sphinx][sphinx] and its extensions can be installed using the
Python package manager via

```console
$ cd <path/to/repository>
# --user installs to a user-specific directory instead of the system
$ pip install --user -r docs/requirements.txt
```

To activate the documentation build targets, the `ACTS_BUILD_DOCS` option has to be set

```console
$ cmake -B <build-dir> -S <path/to/repository> -DACTS_BUILD_DOCS=on
```

Then the documentation can be build with either of the following two build
targets

```console
$ cmake --build <build-dir> docs # default fast option
# or
$ cmake --build <build-dir> docs-with-api # full documentation
```

The default option includes the Doxygen, Sphinx, and the Breathe extension, i.e.
the source code information can be used in the manually written documentation
but the full API documentation is not generated. The second target builds the
full documentation using Exhale to automatically generate the API documentation.
This is equivalent to the public [Read the Docs][rtd_acts] documentation, but
the build takes around ten minutes to finish.

[doxygen]: https://doxygen.nl/
[sphinx]: https://www.sphinx-doc.org
[breathe]: https://breathe.readthedocs.io
[exhale]: https://exhale.readthedocs.io
[rtd_acts]: https://acts.readthedocs.io

## Build options

CMake options can be set by adding `-D<OPTION>=<VALUE>` to the configuration
command. The following command would e.g. enable the unit tests

```console
$ cmake -B <build-dir> -S <source-dir> -DACTS_BUILD_UNITTESTS=ON
```

Multiple options can be given. `cmake` caches the options so that only changed
options must be specified in subsequent calls to configure the project. The
following options are available to configure the project and enable optional
components. q

| Option                                | Description |
|---------------------------------------|-------------|
| ACTS_BUILD_EVERYTHING                 | Build with most options enabled (except HepMC3 and documentation) |
| ACTS_PARAMETER_DEFINITIONS_HEADER     | Use a different (track) parameter definitions header |
| ACTS_BUILD_CUDA_PLUGIN                | Build CUDA plugin |
| ACTS_BUILD_DD4HEP_PLUGIN              | Build DD4hep geometry plugin |
| ACTS_BUILD_DIGITIZATION_PLUGIN        | Build Digitization plugin |
| ACTS_BUILD_IDENTIFICATION_PLUGIN      | Build Identification plugin |
| ACTS_BUILD_JSON_PLUGIN                | Build Json plugin |
| ACTS_USE_BUNDLED_NLOHMANN_JSON        | Use external nlohmann::json installation (on by default) |
| ACTS_BUILD_TGEO_PLUGIN                | Build TGeo plugin |
| ACTS_BUILD_LEGACY                     | Build legacy plugin |
| ACTS_BUILD_FATRAS                     | Build FAst TRAcking Simulation package |
| ACTS_BUILD_EXAMPLES                   | Build standalone examples |
| ACTS_BUILD_EXAMPLES_DD4HEP            | Build DD4hep-based code in the examples |
| ACTS_BUILD_EXAMPLES_GEANT4            | Build Geant4-based code in the examples |
| ACTS_BUILD_EXAMPLES_HEPMC3            | Build HepMC3-based code in the examples |
| ACTS_BUILD_EXAMPLES_PYTHIA8           | Build Pythia8-based code in the examples |
| ACTS_BUILD_BENCHMARKS                 | Build benchmarks |
| ACTS_BUILD_UNITTESTS                  | Build unit tests |
| ACTS_BUILD_INTEGRATIONTESTS           | Build integration tests |
| ACTS_BUILD_DOCS                       | Build documentation |

All Acts-specific options are disabled or empty by default (except for
`ACTS_USE_BUNDLED_NLOHMANN_JSON`) and must be specifically requested. Some of
the options have interdependencies that are automatically handled, e.g. enabling
any of the specific `ACTS_BUILD_EXAMPLES_...` options will also enable the
overall `ACTS_BUILD_EXAMPLES` option. You only need to tell the build system
what you want and it will figure out the rest.

In addition to the Acts-specific options, many generic options are available
that modify various aspects of the build. The following options are some of the
most common ones. For more details, have a look at the annotated list of [useful
CMake variables](https://cmake.org/Wiki/CMake_Useful_Variables) or at the [CMake
documentation](https://cmake.org/documentation/).

| Option               | Description |
|----------------------|-------------|
| CMAKE_BUILD_TYPE     | Build type, e.g. Debug or Release; affects compiler flags |
| CMAKE_CXX_COMPILER   | Which C++ compiler to use, e.g. g++ or clang++ |
| CMAKE_INSTALL_PREFIX | Where to install Acts to |
| CMAKE_PREFIX_PATH    | Search path for external packages |

The build is also affected by some environment variables. They can be set by prepending them to the configuration call:

```console
$ DD4hep_DIR=<path/to/dd4hep> cmake -B <build-dir> -S <source-dir>
```

The following environment variables might be useful.

| Environment variable | Description |
|----------------------|-------------|
| DD4hep_DIR           | Search path for the DD4hep installation |
| HepMC3_DIR           | Search path for the HepMC3 installation |
| Pythia8_DIR          | Search path for the Pythia8 installation |

## Using Acts

When using Acts in your own CMake-based project, you need to include the
following lines in your `CMakeLists.txt` file:

```cmake
find_package (Acts COMPONENTS comp1 comp2 ...)
```

where `compX` are the required components from the Acts project. See the
`cmake` output for more information about which components are available.
