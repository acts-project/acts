## Detector Metadata

The `detray::detectors` component provides the detector metadata headers
(`detray/detectors/<name>_metadata.hpp`). Only hand-written headers, currently
`telescope_metadata.hpp`, live in `include/`. All others are generated at build
time by the `detray_metadata` package in `codegen/`; see
[`codegen/README.md`](codegen/README.md) for the vocabulary, layout and how to
add shipped or external configurations.

The build provides `detray/detectors/metadata_registry.hpp` for every available
metadata family, including telescope and external providers. It defines:

- `DETRAY_FOR_EACH_METADATA(X)`, expanding `X(id)` for each stable family ID;
- `DETRAY_HAS_<ID>_METADATA` feature macros;
- `detray::registry::<id>_metadata<algebra_t>` aliases;
- `detray::registry::metadata_ids` and corresponding `metadata_types<algebra_t>`
  and `detector_types<algebra_t>` tuples, in the same order.

The shipped IDs are `default`, `itk`, `odd`, `toy`, `wire_chamber`, and
`telescope`. IDs describe metadata families, not individual geometries or
algebra/precision choices. The telescope alias fixes the module shape to
`rectangle2D`; other shapes remain usable directly outside this registry.

Availability does not enable backend algorithm instantiations. Traccc still
selects those with `TRACCC_SUPPORTED_DETECTORS`, defaulting to
`default_detector;odd_detector;telescope_detector`. Its aliases, selection
validation, and SYCL tags use this registry.

### External generated metadata

Use an explicit ID to decouple registration from the configuration filename:

```sh
cmake ... '-DDETRAY_GENERATE_METADATA=tracker=/path/to/configuration.py:configure'
```

This generates `detray/detectors/tracker_metadata.hpp` and
`detray::tracker_metadata<algebra_t>`. The older filename-derived selector still
works. The explicit ID is passed to the generator, so the C++ spelling and
metadata name written by IO follow it too.

### External pregenerated metadata

A parent project can register headers supplied by an ordinary CMake target
before adding detray:

```cmake
include(path/to/Detray/cmake/detray-metadata.cmake)
find_package(AtlasMetadata REQUIRED) # provides Atlas::Metadata
# For a locally created target, register on its real target, not its ALIAS.
detray_register_metadata(Atlas::Metadata
    ID tracker
    HEADER atlas/tracker.hpp
    TYPE "atlas::tracker_metadata<@ALGEBRA@, atlas::layout>")
set(DETRAY_METADATA_PROVIDERS Atlas::Metadata)
add_subdirectory(path/to/Detray)
```

The target supplies include paths, dependencies and header installation.
`TYPE` is an algebra-parametric C++ expression; extra template arguments can
be fixed by the provider. Registration properties are exported with locally
built targets. Packages may register before exporting, in which case consumers
only need to list their imported target in `DETRAY_METADATA_PROVIDERS`.
The enclosing project must make external dependency targets available before
loading an installed detray export, just as it does for other CMake dependencies.

IDs must match `[a-z][a-z0-9_]*`. Duplicate IDs, missing provider targets, and
unavailable backend selections fail at configure time. Distinct IDs should
identify distinct C++ metadata types; aliases of the same type are not separate
families.

When ITk moves out of tree, remove its shipped registration and let its provider
register `itk`. Registry consumers and their selections then stay unchanged.
`detray/detectors/io.hpp` supplies `read_registered_detector_json<algebra_t>`,
`read_registered_detector<algebra_t>` (from a payload), and
`write_registered_detector<algebra_t>`. These functions return or accept the
common host detector owner. Writers record the stable family ID in every JSON
component. Readers also accept the older `<id>_metadata` spelling. Missing or
unknown IDs throw; pass an explicit ID to read legacy data. Detector display
names do not select metadata types.

The configure/export checks run with `ctest -R detray_metadata_registry` in a
build with `DETRAY_BUILD_TESTING=ON`.

### Host type erasure

`detray/core/host_detector.hpp` and `detray/core/detector_buffer.hpp` provide
move-only owners for concrete detectors and buffers. They belong to `detray::core`
and do not depend on the registry: external concrete types can also be stored.
The memory resource must outlive the owner. Buffer views borrow the buffer;
asynchronous work must finish before the buffer or its resource is destroyed.

Operations explicitly choose a tuple of supported concrete detector types:

```cpp
detray::host_detector detector;
detector.set<my_metadata>(std::move(concrete_detector));
auto count = detray::host_detector_visitor<my_detector_types>(
    detector, []<typename T>(const T& concrete) {
        return concrete.volumes().size();
    });
auto buffer = detray::buffer_from_host_detector<my_detector_types>(
    detector, resource, copy);
```

`detector_buffer_visitor` similarly calls the operation with a concrete view.
Dispatch happens on the host; kernels and navigation still use concrete types.
An empty owner or a type outside the operation's supported list throws
`std::invalid_argument`. Visitor return types must agree across the list;
references are preserved and borrow their underlying owner.
