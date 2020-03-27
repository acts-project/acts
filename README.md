# A Common Tracking Software (Acts) Project

[![coverage](https://badgen.net/codecov/c/github/acts-project/acts/master)](https://codecov.io/gh/acts-project/acts/branch/master) [![10.5281/zenodo.3606011](https://badgen.net/badge/DOI/10.5281%2Fzenodo.3606011/blue)](https://doi.org/10.5281/zenodo.3606011) [![Chat on Mattermost](https://badgen.net/badge/chat/on%20mattermost/cyan)](https://mattermost.web.cern.ch/acts/)

1. [Introduction](#intro)
2. [Mailing list](#mailing-list)
6. [License and authors](#license-authors)

# <a name="intro">Introduction</a>

This project is an experiment-independent set of track reconstruction tools. The
main philosophy is to provide high-level track reconstruction modules that can
be used for any tracking detector. The description of the tracking detector's
geometry is optimized for efficient navigation and quick extrapolation of
tracks. Converters for several common geometry description languages exist.
Having a fast, yet largely customizable implementation of track
reconstruction algorithms was a primary objective for the design of this
toolset. Additionally, the applicability to real-life HEP experiments played a
major role in the development process. Apart from algorithmic code, this project
also provides an event data model for the description of track parameters and
measurements.

Key features of this project include:

*   Tracking geometry description which can be constructed from TGeo, DD4Hep,
    or customized input.
*   Simple and efficient event data model.
*   Fast and highly flexible algorithms for track propagation and fitting.
*   Seed finding algorithms.
*   Fast track simulation tools based on the reconstruction geometry
    description.

The git repository for the Acts project can be found at <a href="https://gitlab.cern.ch/acts/acts-core.git">https://gitlab.cern.ch/acts/acts-core.git</a>.

To get started, you can have a look at our [Getting Started Guide](getting_started.md)

# <a name="mailing-list">Mailing list</a>

In order to receive the latest updates, users of the Acts project are encouraged to subscribe to [acts-users@cern.ch](https://e-groups.cern.ch/e-groups/Egroup.do?egroupName=acts-users). This list provides:
- regular updates on the software,
- access to the Acts JIRA project for bug fixes/feature requests,
- a common place for asking any kind of questions.

If you find a bug, have a feature request, or want to contribute to Acts, please have a look at the [contribution guide](CONTRIBUTING.md).

# <a name="license-authors">License and authors</a>

This project is published under the terms of the Mozilla Public License, v. 2.0.
A copy of the license can be found in the [LICENSE](LICENSE) file or at
http://mozilla.org/MPL/2.0/ .

Contributors to the Acts project are listed in the [AUTHORS](AUTHORS.md) file.

The Acts project is based on the ATLAS tracking software. A list of contributors
to the ATLAS tracking repository can be found at
http://acts.web.cern.ch/ACTS/ATLAS_authors.html .

The Acts project contains copies of the following external packages:

*   [dfelibs](https:://gitlab.cern.ch/msmk/dfelibs) by Moritz Kiehn licensed
    under the MIT license.
*   [JSON for Modern C++](https://github.com/nlohmann/json) by Niels Lohmann
    licensed under the MIT License.
