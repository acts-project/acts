# Detray Python Package

## Detector Module

This module contains the tools to generate the detector metadata headers used to define a detray detector type. This can be steered manually or automatically during the cmake configuration step.

### Generating custom Detector Metadata

The detector metadata is a C++ `struct` used by the detray core library as a template parameter to gather the compile-time information on the capabilities of the detector modelling, for example:

`using odd_detector_t = detray::host::detector<detray::odd_metadata<algebra_t>>`,

with the `algebra_t` defined according to the linear algebra backend to be used, for instance the `detray::array<float>` plugin.

A generic metadata type (the default metadata) can represent all supported
detectors, at the expense of increased build times and potentially slower
reconstruction. Detector-specific configurations live in
`detray/detectors/configurations` within this Python package.

Select metadata headers through CMake, for example:

```shell
cmake -S detray -B detray-build -DDETRAY_GENERATE_METADATA=odd_metadata
```

The codegen manifest invokes the common generator with the selected detector
configuration, or supplies a prebuilt header. To generate a header manually:

```shell
source detray-build/python/setup.sh
python3 detray/detectors/python/generate_metadata.py odd_metadata.hpp --detector odd.add_odd_types
```

### Metadata Configuration Manual

Generate a metadata representation and dump the header file:

```python
from detray.detectors import metadata, metadata_generator

md = metadata("detector_name")

[...]

metadata_generator(md)
```

Adding surface shapes for passive or sensitive detector elements or portals between detector sub-volumes:
