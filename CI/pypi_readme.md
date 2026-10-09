[![ACTS logo](https://raw.githubusercontent.com/acts-project/acts/main/docs/figures/acts_logo_colored.svg)](https://github.com/acts-project/acts)
# ACTS Common Tracking Software

or *A Common Tracking Software* if you do not like recursive acronyms

[![10.5281/zenodo.5141418](https://zenodo.org/badge/DOI/10.5281/zenodo.5141418.svg)](https://doi.org/10.5281/zenodo.5141418)
[![Chat on Mattermost](https://badgen.net/badge/chat/on%20mattermost/cyan)](https://mattermost.web.cern.ch/acts/)
[![Latest release](https://badgen.net/github/release/acts-project/acts)](https://github.com/acts-project/acts/releases)

ACTS is an experiment-independent toolkit for (charged) particle track
reconstruction in (high energy) physics experiments implemented in modern C++.

More information can be found in the [ACTS documentation](https://acts-project.github.io/).

## Installation

The easiest way to get the ACTS Python bindings is via [PyPI](https://pypi.org/project/pyacts/):

```console
pip install pyacts
```

The `pyacts` distribution includes bindings for the core library, Fatras, and
the Examples framework.

> **Warning**
> The PyPI package does not include all optional plugins (e.g. DD4hep/ODD, Geant4, ROOT).
> To use those, build from source with the relevant CMake options enabled.

> **Note**
> Even though the name of the PyPI package is `pyacts`, the python module must be imported
> as `import acts`.


The wheel includes `.pyi` type stubs and supports Linux and macOS with Python 3.11 or newer.

For a **source build**, enable the Python bindings and Examples framework, then load the generated
`<build>/this_acts_withdeps.sh` script before importing `acts`:

```console
cmake -B <build> -S <source> -DACTS_BUILD_PYTHON_BINDINGS=ON -DACTS_BUILD_EXAMPLES=ON
cmake --build <build> --target ActsPythonBindings
source <build>/this_acts_withdeps.sh
```

Enable optional components with their CMake options. For ROOT readers, writers, and evaluation,
use `-DACTS_BUILD_EXAMPLES_ROOT=ON`; this also enables the ROOT plugin. A source build needs Python
development headers. The `ActsPythonBindings` target builds the enabled bindings, and the setup
script adds the built modules and dependencies to your environment.

| Capability | PyPI (`pyacts`) | Full installation (source build) |
| :--- | :--- | :--- |
| Plugins | Arrow, JSON | All available plugins, including ROOT, DD4hep, or Geant4 |
| I/O | CSV, Parquet, particle and sim-hit ROOT `uproot` readers | Additionally native ROOT and EDM4hep |
| Geometry | `GenericDetector`, Tracking geometry from JSON | Directly build detector geometries with DD4hep or TGeo |
| Evaluation | In-memory performance writers and plotting in Python | Additionally native ROOT performance writers |

## Getting started

An `acts.examples.Sequencer` runs readers, algorithms, and writers over events. This first script
generates muons and prints them; it works with the PyPI wheel:

```python
import acts
from acts import UnitConstants as u
import acts.examples
from acts.examples.simulation import (
    EtaConfig,
    MomentumConfig,
    ParticleConfig,
    addParticleGun,
)

s = acts.examples.Sequencer(events=5, numThreads=1)
addParticleGun(
    s,
    momentumConfig=MomentumConfig(1 * u.GeV, 10 * u.GeV, transverse=True),
    etaConfig=EtaConfig(-2.0, 2.0),
    particleConfig=ParticleConfig(1, acts.PdgParticle.eMuon),
    rnd=acts.examples.RandomNumbers(seed=42),
    printParticles=True,
)
s.run()
```

For convenience, the package contains helpers in `acts.examples.simulation` and `acts.examples.reconstruction` that add entire workflow steps with `addXYZ`-functions (such as simulation, seeding, ...)
and connect their data keys. For full flexibility, the underlying sequencer algorithms can be directly configured as well. For simple but complete chains, see the
[PyPI finding and fitting demo](https://github.com/acts-project/acts/blob/main/Examples/Scripts/Python/pypi_finding_fitting_demo.py)
or the [truth-tracking Kalman example](https://github.com/acts-project/acts/blob/main/Examples/Scripts/Python/truth_tracking_kalman.py)
(source build). Their helper signatures and default keys are in
[simulation.py](https://github.com/acts-project/acts/blob/main/Python/Examples/python/simulation.py)
and [reconstruction.py](https://github.com/acts-project/acts/blob/main/Python/Examples/python/reconstruction.py).
