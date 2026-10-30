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

### External metadata configurations

A configuration outside this repository can provide a function accepting a
`metadata` object, just like the packaged configurations:

```python
# /path/to/custom.py
from detray.detectors.configurations.toy import add_toy_types

def configure(md):
    add_toy_types(md)
    # Add or replace types here.
```

Select it with
`-DDETRAY_GENERATE_METADATA='/path/to/custom.py:configure'` (semicolon-separated
with any shipped metadata names), or invoke `generate_metadata.py` directly with
`--detector /path/to/custom.py:configure`. The file stem names the detector and
its installed header, `detray/detectors/custom_metadata.hpp`; it must be a Python
identifier and must not reuse a shipped metadata name. CMake tracks changes to
the configuration file. Imports should come from the detray package or other
packages available in the generator environment.

External configurations use the common generator runtime and prebuilt lookup,
but are not included in release pregeneration. They generate locally unless a
header is supplied at
`<prebuilt-dir>/Detray/detectors/detray/detectors/custom_metadata.hpp`.
Strict prebuilt mode rejects a missing header. Supplying and keeping that custom
prebuilt header up to date is the downstream project's responsibility.
