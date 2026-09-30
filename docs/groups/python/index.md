@defgroup python_bindings Python Bindings
@brief Use ACTS from Python, with the PyPI wheel or a source build.

<span id="python-bindings-overview"></span>The `pyacts` distribution is imported as `acts`.
It exposes the core library, Fatras, and the
Examples framework. A source build can also expose optional plugins.

## Installation

Install the [PyPI package](https://pypi.org/project/pyacts/) with:

```console
pip install pyacts
```

The wheel includes `.pyi` type stubs and supports Linux and macOS with Python 3.11 or newer.

For a source build, enable the Python bindings and Examples framework, then load the generated
`<build>/this_acts_withdeps.sh` script before importing `acts`:

```console
cmake -B <build> -S <source> -DACTS_BUILD_PYTHON_BINDINGS=ON -DACTS_BUILD_EXAMPLES=ON
cmake --build <build>
source <build>/this_acts_withdeps.sh
```

Enable optional plugins with their CMake options (for example,
`-DACTS_BUILD_PLUGIN_ROOT=ON`). A source build needs Python development headers.

| Capability | PyPI (`pyacts`) | Full installation (source build) |
| :--- | :--- | :--- |
| Plugins | Arrow, JSON, and FPE monitoring; no ROOT, DD4hep, or Geant4 | Optional plugins enabled through CMake |
| I/O | CSV, Parquet, and HepMC3; ROOT files via Python `uproot` readers | Same formats, plus ROOT and EDM4hep when enabled |
| Geometry | Load tracking geometry from JSON; GenericDetector examples | Can also build detector geometry with DD4hep/ODD when enabled |
| Evaluation | In-memory performance writers and plotting in Python | Same Python writers, plus ROOT performance writers when enabled |

Optional features in the full installation depend on which components and dependencies you build.

## Getting started

An `acts.examples.Sequencer` runs readers, algorithms, and writers over events. This first script
generates muons and prints them; it works with the PyPI wheel:

@snippet{trimleft} examples/test_python_getting_started.py First Python run

The helpers in `acts.examples.simulation` and `acts.examples.reconstruction` add later stages
and connect their data keys. For complete chains, see the
[PyPI finding and fitting demo](https://github.com/acts-project/acts/blob/main/Examples/Scripts/Python/pypi_finding_fitting_demo.py)
or the [truth-tracking Kalman example](https://github.com/acts-project/acts/blob/main/Examples/Scripts/Python/truth_tracking_kalman.py)
(source build). Their helper signatures and default keys are in
[simulation.py](https://github.com/acts-project/acts/blob/main/Python/Examples/python/simulation.py)
and [reconstruction.py](https://github.com/acts-project/acts/blob/main/Python/Examples/python/reconstruction.py).

## Next steps

- @ref python_custom_algorithms "Custom algorithms" — add Python readers and algorithms.
- @ref python_data_io "Data reading and writing" — choose event and geometry formats.
- @ref python_performance_plotting "Performance evaluation and plotting" — inspect and plot results.
