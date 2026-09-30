@defgroup python_bindings Python Bindings
@brief Use ACTS from Python, with the PyPI wheel or a source build.

The `pyacts` distribution is imported as `acts`. It exposes the core library, Fatras, and the
Examples framework. A source build can also expose optional plugins.

| Capability | PyPI (`pyacts`) | Full installation (source build) |
| :--- | :--- | :--- |
| Plugins | Arrow, JSON, and FPE monitoring; no ROOT, DD4hep, or Geant4 | Optional plugins enabled through CMake |
| I/O | CSV and Parquet; ROOT files via the Python `uproot` readers | CSV and Parquet, plus ROOT, EDM4hep, or HepMC3 when enabled |
| Geometry | Load tracking geometry from JSON; GenericDetector examples | Can also build detector geometry with DD4hep/ODD when enabled |
| Evaluation | In-memory performance writers and plotting in Python | Same Python writers, plus ROOT performance writers when enabled |

The wheel includes `.pyi` type stubs and supports Linux and macOS with Python 3.11 or newer.
Optional features in the full installation depend on which components and dependencies you build.

## Installation

Install the [PyPI package](https://pypi.org/project/pyacts/) with:

```console
pip install pyacts
```

For a source build, enable the Python bindings and Examples framework, then load the generated
`<build>/this_acts_withdeps.sh` script before importing `acts`:

```console
cmake -B <build> -S <source> -DACTS_BUILD_PYTHON_BINDINGS=ON -DACTS_BUILD_EXAMPLES=ON
cmake --build <build>
source <build>/this_acts_withdeps.sh
```

Enable optional plugins with their CMake options (for example,
`-DACTS_BUILD_PLUGIN_ROOT=ON`). A source build needs Python development headers.

## Getting started

An `acts.examples.Sequencer` runs readers, algorithms, and writers over events. The helpers in
`acts.examples.simulation` and `acts.examples.reconstruction` add complete stages and connect
their data keys. A typical script creates a sequencer, adds generation, simulation, digitization,
and reconstruction, then calls `s.run()`.

Start with `Examples/Scripts/Python/pypi_finding_fitting_demo.py` for a wheel-compatible chain,
or `Examples/Scripts/Python/truth_tracking_kalman.py` for a source-build example. The helper
signatures and default data keys are in `Python/Examples/python/simulation.py` and
`reconstruction.py`.

## Next steps

- @ref python_custom_algorithms "Custom algorithms" — add Python readers and algorithms.
- @ref python_data_io "Data reading and writing" — choose event and geometry formats.
- @ref python_performance_plotting "Performance evaluation and plotting" — inspect and plot results.
