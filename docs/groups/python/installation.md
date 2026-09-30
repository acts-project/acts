@defgroup python_installation Installation
@ingroup python_bindings
@brief Getting the `acts` Python package, from PyPI or from source.

# Installation from PyPI

The easiest way to get the ACTS Python bindings is via [PyPI](https://pypi.org/project/pyacts/):

```console
pip install pyacts
```

A nightly development build is also available from [TestPyPI](https://test.pypi.org/project/pyacts/):

```console
pip install -i https://test.pypi.org/simple pyacts
```

> [!note]
> Even though the name of the PyPI package is `pyacts`, the python module must be imported
> as `import acts`.

## What is in the wheel

The `pyacts` distribution includes bindings for the core library, Fatras, and the Examples
framework, built with a fixed configuration (the `python-wheel` CMake preset):

- **Arrow/Parquet I/O is built in** (`acts.examples.arrow`, `ParquetReader`/`ParquetWriter`) —
  this is the wheel's primary route for reading and writing event data; see
  @ref python_data_io.
- **ROOT is not built in.** `acts.root`, `acts.examples.root`, and every `Root*Reader`/`Root*Writer`
  are unavailable. Use the ROOT-free readers (`acts.examples.uproot`), the Arrow/Parquet path, or
  the ROOT-free performance writers instead (@ref python_data_io, @ref python_performance_plotting).
- **DD4hep, the OpenDataDetector, and Geant4 are not built in.** Load geometry from JSON
  instead (`acts.json.TrackingGeometryJsonConverter`), or build from source if you need DD4hep.
- **Type stubs (`.pyi`) are included**, generated with `pybind11-stubgen`, so IDEs get
  autocompletion without a source build.
- Supported on Linux (manylinux_2_34) and macOS, Python ≥ 3.11.

> [!warning]
> The PyPI package does not include DD4hep/ODD, Geant4, or ROOT. To use those, build from
> source with the relevant CMake options enabled.

`Examples/Scripts/Python/pypi_finding_fitting_demo.py` is a self-contained demo that only uses
what the wheel provides — a good starting point to check your installation, and the reference
example for @ref python_custom_algorithms.

# Building from source

To use the full Python bindings including optional components such as DD4hep/ODD, Geant4, ROOT,
build ACTS from source with `ACTS_BUILD_PYTHON_BINDINGS=ON`.

```console
cmake -B <build> -S <source> -DACTS_BUILD_PYTHON_BINDINGS=ON
```

Bindings for additional components (plugins, Examples framework) are built automatically if
they are enabled.

```console
cmake -B <build> -S <source> -DACTS_BUILD_PYTHON_BINDINGS=ON -DACTS_BUILD_EXAMPLES=ON
```

Building requires a Python installation including the development headers.
You can then build the special target `ActsPythonBindings` to build everything
that can be accessed from Python:

```console
cmake --build <build> -- ActsPythonBindings
```

The build creates a setup script `$BUILD/this_acts_withdeps.sh` which modifies
`$PYTHONPATH` so that you can import the `acts` module in Python.

> [!warning]
> The old flag `ACTS_BUILD_EXAMPLES_PYTHON_BINDINGS` is **deprecated**.
> Use `ACTS_BUILD_PYTHON_BINDINGS` (together with `ACTS_BUILD_EXAMPLES` if
> needed) instead.

To generate the `.pyi` type stubs in a source build as well, add `-DACTS_GENERATE_PYTHON_STUBS=ON`
(requires [uv](https://docs.astral.sh/uv/)).
