# Detray Metadata Generator

`detray_metadata` generates the detector metadata headers that define a detray
detector type. It is a build-time code generator, independent of the `detray`
Python bindings in `Detray/python`.

## Vocabulary

| Term | Meaning | Example |
|---|---|---|
| name | the key used in CMake and on the command line | `toy` |
| configuration | Python module `<name>_detector` with a `configure(md)` function | `detray_metadata/configurations/toy_detector.py` |
| metadata | the generated C++ type list | `detray::toy_metadata<algebra_t>` in `detray/detectors/toy_metadata.hpp` |
| detector type | metadata combined with an algebra | `detray::host::detector<detray::toy_metadata<algebra_t>>`, traccc's `toy_detector` |

The detector metadata is a C++ `struct` used by the detray core library as a
template parameter to gather the compile-time information on the capabilities
of the detector modelling, with the `algebra_t` defined according to the linear
algebra backend to be used, for instance the `detray::array<float>` plugin.

A generic metadata type (`default`) can represent all supported detectors, at
the expense of increased build times and potentially slower reconstruction.

## Layout

```
Detray/detectors/
├── include/detray/detectors/      # hand-written metadata headers only (telescope)
├── CMakeLists.txt                 # generates the shipped and external headers
└── codegen/
    ├── generate_metadata.py       # CLI used by the codegen manifest
    └── detray_metadata/
        ├── impl/                  # metadata model and header rendering
        ├── options/               # command line options
        ├── utils/                 # reusable tracker/material defaults
        └── configurations/        # <name>_detector.py, one per shipped detector
```

## Building

The shipped configurations (`default`, `odd`, `toy`, `wire_chamber`) are always
generated through the ACTS codegen manifest, or taken from prebuilt headers.
Generated headers are never committed.

The ITk metadata is not shipped: traccc generates it as an external
configuration from `Traccc/metadata/itk_detector.py`, selected by
`TRACCC_ITK_METADATA_CONFIGURATION`. It is therefore not part of release
pregeneration. The configuration filename must be `itk_detector.py` or `itk.py`
so that it generates the `itk_metadata.hpp` header expected by traccc. For a
release archive build with `ACTS_CODEGEN_REQUIRE_PREBUILT=ON`, supply
`Detray/detectors/detray/detectors/itk_metadata.hpp` under
`ACTS_CODEGEN_PREBUILT_DIR`, alongside the shipped prebuilt headers. Otherwise,
turn strict mode off to generate the ITk header locally.

To generate a header manually:

```shell
PYTHONPATH=Detray/detectors/codegen \
    python3 Detray/detectors/codegen/generate_metadata.py odd_metadata.hpp --detector odd
```

### Adding a shipped detector

Add `configurations/<name>_detector.py` with a `configure(md)` function, a
manifest unit `Detray/detectors/detray/detectors/<name>_metadata.hpp` with
`"arguments": ["--detector", "<name>"]`, and `<name>` to the list in
`Detray/detectors/CMakeLists.txt`.

### External configurations

A configuration outside this repository follows the same convention:

```python
# /path/to/custom_detector.py
from detray_metadata import add_silicon_tracker_defaults


def configure(md):
    add_silicon_tracker_defaults(
        md, use_homogeneous_mat=True, use_mat_maps=True, add_trapezoid=True
    )
```

Select it with
`-DDETRAY_GENERATE_METADATA='/path/to/custom_detector.py:configure'`
(semicolon-separated), or pass `--detector /path/to/custom_detector.py:configure`
to `generate_metadata.py`. A trailing `_detector` is dropped from the file stem
to form the name, so this produces `detray/detectors/custom_metadata.hpp` with
`detray::custom_metadata`. The name must not reuse a shipped one. CMake tracks
changes to the configuration file. Each name can be generated only once per
build, so do not also pass the ITk configuration via `DETRAY_GENERATE_METADATA`
when building traccc.

External configurations use the common generator runtime and prebuilt lookup,
but are not part of release pregeneration. They generate locally unless a header
is supplied at
`<prebuilt-dir>/Detray/detectors/detray/detectors/custom_metadata.hpp`; strict
prebuilt mode rejects a missing header.

Alternatively, `detray_metadata` can be imported from a script of your own that
builds a `metadata` object and calls `metadata_generator` directly.

### Metadata Configuration Manual

Generate a metadata representation and dump the header file:

```python
from detray_metadata import metadata, metadata_generator

md = metadata("detector_name")

[...]

metadata_generator(md)
```

Adding surface shapes for passive or sensitive detector elements or portals between detector sub-volumes:
