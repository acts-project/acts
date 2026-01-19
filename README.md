# ACTS Common Tracking Software

or *A Common Tracking Software* if you do not like recursive acronyms

[![10.5281/zenodo.5141418](https://zenodo.org/badge/DOI/10.5281/zenodo.5141418.svg)](https://doi.org/10.5281/zenodo.5141418)
[![Chat on Mattermost](https://badgen.net/badge/chat/on%20mattermost/cyan)](https://mattermost.web.cern.ch/acts/)
[![Coverage](https://sonarcloud.io/api/project_badges/measure?project=acts-project_acts&metric=coverage)](https://sonarcloud.io/summary/new_code?id=acts-project_acts)
[![Latest release](https://badgen.net/github/release/acts-project/acts)](https://github.com/acts-project/acts/releases)
[![Status](https://badgen.net/github/checks/acts-project/acts/main)](https://github.com/acts-project/acts/actions)
[![Metrics](https://badgen.net/badge/metric/tracker/purple)](https://acts-project.github.io/metrics/)

ACTS is an experiment-independent toolkit for (charged) particle track
reconstruction in (high energy) physics experiments implemented in modern C++.

More information can be found in the [ACTS documentation](https://acts.readthedocs.io/).

## Quick start

ACTS is developed in C++ and is build using [CMake](https://cmake.org). Building
the core library requires a C++20 compatible compiler,
[Boost](https://www.boost.org), and [Eigen](https://eigen.tuxfamily.org). The
following commands will clone the repository, configure, and build the core
library

```sh
git clone https://github.com/acts-project/acts <source>
cmake -B <build> -S <source>
cmake --build <build>
```

For more details, e.g. specific versions and additional dependencies, have a
look at the [getting started guide](docs/getting_started.md). If you find a bug,
have a feature request, or want to contribute to ACTS, have a look at the
[contribution guidelines](CONTRIBUTING.md).

## Versioning and public API

Release versions follow [semantic versioning](https://semver.org/spec/v2.0.0.html)
to indicate whether a new version contains breaking changes within the public API.
Currently, only a limited part of the visible API is considered the public API
and subject to the semantic versioning rules. The details are outlined in the
[versioning and public API documentation](docs/versioning.rst).

## Repository organization

The repository contains all code of the ACTS projects, not just the core library
that a physics experiment is expected to use as part of its reconstruction code.
All optional components are disabled by default. Please see the
[getting started guide](docs/getting_started.md) on how-to enable them.


| *Folder*     | *include path*                              | *Namespace*       | *Python module*          |
| :---         | :---                                        | :---              | :---                     |
| `Core/`      | `#include "Acts/..."`                       | `Acts::`          | `acts`                   |
| `Plugins/`   | `#include "ActsPlugins/PluginName/..."`     | `ActsPlugins::`   | `acts.pluginname`        |
| `Fatras/`    | `#include "ActsFatras/..."`                 | `ActsFatras::`    | `acts.fatras`            |
| `Alignment/` | `#include "ActsAlignment/..."`              | `ActsAlignment::` | `acts.alignment`         |
| `Examples/`  | `#include "ActsExamples/..."`               | `ActsExamples::`  | `acts.examples`          |
| `Python/`    | `#include "ActsPython/..."` (not much used) | `ActsPython::`    |  N/A                     |
| `Tests/`     | `#include "ActsTests/..."` (not much used)  | `Acts::Tests::`   |  N/A                     |


Short summary of the modules and directories:
-   `Core/` contains all the core functionality with minimal dependencies and no framework or I/O related code
-   `Plugins/` contains plugins for core functionality that require
    additional external packages.
-   `Fatras/` provides fast track simulation tools based on the core
    library.
-   `Examples/` contains simulation and reconstruction examples. These are
    internal tools for manual full-chain development and tests.
-   `Tests/` contains automated unit tests, integration tests, and
    (micro-)benchmarks.
-   `Python/` contains the python bindings for the different modules
-   `thirdparty/` contains external dependencies that are usually not available
    through the system package manager.

## Authors and license

Contributors to the ACTS project are listed in the [AUTHORS](AUTHORS) file.

The ACTS project is published under the terms of the Mozilla Public License, v. 2.0.
A copy of the license can be found in the [LICENSE](LICENSE) file or at
https://mozilla.org/MPL/2.0/ .

The ACTS project contains copies of the following external packages:

-   [OpenDataDetector](https://github.com/acts-project/OpenDataDetector)
    licensed under the MPLv2 license.
