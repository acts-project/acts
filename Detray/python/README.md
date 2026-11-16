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

## Detector API

`detray.Detector` is the common Python type for all registered metadata families.
`detray.core.availableMetadata` lists the families available in this build.
There is one configured scalar/algebra choice per build.

```python
import detray
import detray.io

cfg = detray.io.DetectorReaderConfig().addFile("detector_geometry.json")
detector, names = detray.io.readDetector(detray.HostMemoryResource(), cfg)
print(detector.metadataId, detector.name(names))
print(len(detector.volumes), len(detector.surfaces))
print(detector.surfaces[0].identifier)
```

The reader selects the family from the file header. Use `metadataId="odd"`
(or another registered ID) explicitly for older files without metadata.
Unknown IDs raise `ValueError`. `writeDetector` accepts the same common type
and records its stable metadata ID.

`volumes`, `surfaces`, and `portals` are read-only sequences with negative-index
support. They retain the detector and its allocator. Their `VolumeInfo` and
`SurfaceInfo` elements are independent snapshots containing indices, surface
ranges, identifiers and surface classifications. Metadata-specific mask,
material and accelerator stores are not part of this initial inspection API.
`VolumeGraph` captures a printable/DOT snapshot of the detector connectivity.

Detectors returned by IO and test builders keep their memory resource alive.
SVG illustrators retain their detector and names. Standalone detray imports
never import ACTS; ACTS bridge modules use these canonical registrations.

## ACTS bridge and propagation

`acts.detray.DetrayGeometryConverter.convert(resource, context, geometry,
metadataId="odd")` returns a result whose `detector` is this same
`detray.Detector`. The extracted detector keeps the conversion result, allocator,
and ACTS tracking geometry alive. `acts.vecmem` forwards the detray resource types.

Use `detray.io.readDetector` and `writeDetector` for both standalone and converted
geometries. `acts.examples.detray.StraightLinePropagator` accepts the common
`Detector` and dispatches to its registered concrete type; it retains the detector
and propagation allocator. The validation scripts accept `--detray-metadata` for
explicit metadata selection, including legacy JSON files.

The ACTS bridge tests in `Python/Plugins/tests/test_detray.py` cover import order,
conversion, JSON round trips and allocator ownership. ACTS jobs discover these
alongside detray's own tests.
