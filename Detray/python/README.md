# Python bindings

The `detray` package can be built through either the ACTS top-level project or
`cmake -S Detray` in the same checkout. It does not require importing or building
ACTS. Both entry points use `cmake/ActsPython.cmake` for extension creation,
installation, relative library paths, Python file staging, and environment setup,
and `cmake/ActsPybindStubs.cmake` for stub generation. Package-specific CMake files
supply the sources, dependencies, destinations, and feature settings.

For a minimal standalone build:

```sh
cmake -S Detray -B build-detray -DDETRAY_BUILD_PYTHON_BINDINGS=ON
cmake --build build-detray --target DetrayPythonBindings DetrayIoPythonBindings
source build-detray/python/setup.sh
python -m pytest Detray/python/Core/tests Detray/python/Io/tests/test_io_configs.py
```

Use the Python interpreter selected by CMake and install `pytest` and `numpy` in
that environment. `DETRAY_GENERATE_PYTHON_STUBS=ON` enables stub generation. The
standalone install destination defaults to `python/detray` under the prefix and
can be overridden with `DETRAY_PYTHON_INSTALL_DIR`.

In the ACTS build, Python bindings plus the traccc plugin enable detray bindings;
the ACTS stub-generation setting applies to detray as well. The wheel layout
installs the top-level `acts` and `detray` packages alongside each other.

The full test suite additionally uses SVG tools and validation utilities:

```sh
cmake --preset detray-full-fp32 -S Detray -B build-detray-full \
    -DDETRAY_BUILD_PYTHON_BINDINGS=ON
cmake --build build-detray-full
source build-detray-full/python/setup.sh
python -m pytest Detray/python
```

Core and I/O configuration tests deliberately avoid those optional utilities.
The standalone CI jobs run the full suite. ACTS Python jobs run the same baseline
tests together with ACTS core basics in a separate pytest invocation, covering
source-path shadowing as well as imports. Missing required modules fail collection;
optional suites are selected explicitly rather than hiding import failures.

ACTS interoperability tests belong under `Python`; tests that only exercise
detray belong here. Shared helper changes must trigger both build entry points.
