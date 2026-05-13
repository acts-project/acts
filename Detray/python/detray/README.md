# Detray Python Package

## Detector Module

This module contains the tools to generate the detector metadata headers used to define a detray detector type. This can be steered manually or automatically during the cmake configuration step.

### Generating custom Detector Metadata

The detector metadata is a C++ `struct` used by the detray core library as a template parameter to gather the compile-time information on the capabilities of the detector modelling, for example:

`using odd_detector_t = detray::detector<detray::odd_metadata<algebra_t>>`,

with the `algebra_t` defined according to the linear algebra backend to be used, for instance the `detray::array<float>` plugin.

A generic metadata type (the default metadata) exists, which comprises most of the geometry modelling capabilities detray has to offer, however, at the expense of increased build times and likely higher runtime of subsequent track reconstruction pipelines. A custom metadata can be defined by providing a python script that configures the detray metadata code generator,
which can be invoked during the configuration stage of the cmake build:
```shell
cmake -S detray -B detray-build -DDETRAY_METADATA_GENERATOR=path/to/custom_metadata.py
```
The resulting metadata header, ready to be used in downstream projects, will be available in the `detray/detectors` folder. Example metadata generation scripts, for instance for the ACTS Open Data Detector (ODD detector mentioned above), can be found in the `detray/detectors/python/` folder. They can be invoked manually by:
```shell
source detray-build/python/setup.sh
python3 detray/detectors/python/odd_metadata.py
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
