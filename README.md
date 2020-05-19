# Acts Common Tracking Software

or *A Common Tracking Software* if you do not like recursive acronyms

[![10.5281/zenodo.3606011](https://badgen.net/badge/DOI/10.5281%2Fzenodo.3606011/blue)](https://doi.org/10.5281/zenodo.3606011)
[![Chat on Mattermost](https://badgen.net/badge/chat/on%20mattermost/cyan)](https://mattermost.web.cern.ch/acts/)
[![coverage](https://badgen.net/codecov/c/github/acts-project/acts/master)](https://codecov.io/gh/acts-project/acts/branch/master)

Acts is an experiment-independent toolkit for (charged) particle track
reconstruction in (high energy) physics experiments implemented in modern C++.

More information can be found on the [Acts website](https://cern.ch/acts) and
in the [Acts documentation](https://acts.readthedocs.io/).

## Using Acts

Acts is developed in C++ and is build using [CMake](https://cmake.org). Building
the core library requires a C++17 compatible compiler,
[Boost](http://boost.org), and [Eigen](http://eigen.tuxfamily.org). The
following commands will clone the repository, configure, and build the core
library

```sh
git clone https://github.com/acts-project/acts <source-dir>
cmake -B <build-dir> -S <source-dir>
cmake --build <build-dir>
```

For more details, e.g. specific versions and additional dependencies, have a
look at the [getting started guide](docs/getting_started.md). If you find a bug,
have a feature request, or want to contribute to Acts, have a look at the
[contribution guide](CONTRIBUTING.md).

## Repository organization

The repository contains all code of the Acts projects, not just the core library
that a physics experiment is expected to use as part of its reconstruction code.
All optional components are disabled by default. Please see the
[getting started guide](docs/getting_started.md) on how-to enable them.

-   `Core/` contains the core library that provides functionality in the `Acts`
    namespace.
-   `Plugins/` contains plugins for core functionality that requires
    additional external packages. The functionality also resides in the `Acts`
    namespace.
-   `Fatras/` provides fast track simulation tools based on the core
    library. This is not part of the core functionality and thus resides in the
    separate `ActsFatras` namespace.
-   `Examples/` contains simulation and reconstruction examples. These are
    internal tools for manual full-chain development and tests and reside in
    the `ActsExamples` namespace.
-   `Tests/` contains automated unit tests, integration tests, and
    (micro-)benchmarks.
-   `thirdparty/` contains external dependencies that are usually not available
    through the system package manager.

## Authors and license

Contributors to the Acts project are listed in the [AUTHORS](AUTHORS) file.

The Acts project is published under the terms of the Mozilla Public License, v. 2.0.
A copy of the license can be found in the [LICENSE](LICENSE) file or at
http://mozilla.org/MPL/2.0/ .

The Acts project contains copies of the following external packages:

-   [dfelibs](https:://gitlab.cern.ch/msmk/dfelibs) by Moritz Kiehn
    licensed under the MIT license.
-   [JSON for Modern C++](https://github.com/nlohmann/json) by Niels Lohmann
    licensed under the MIT License.
