## Generating Detector Metadata

`python/generate_metadata.py` generates metadata for the detector configurations
in `../python/detray/detectors/configurations`. Select the required headers during
CMake configuration, for example `-DDETRAY_GENERATE_METADATA="odd_metadata"`.
The common codegen manifest supplies the detector argument and supports prebuilt
headers.

To generate a header manually:

```shell
source path-to-detray-build/python/setup.sh
python3 detray/detectors/python/generate_metadata.py odd_metadata.hpp --detector odd.add_odd_types
```

To add a detector, define its configuration function in the package and reference
it as `module.function` in the codegen manifest. The module name determines the
detector name; no generator registry needs updating.
