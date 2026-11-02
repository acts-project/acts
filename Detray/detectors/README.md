## Generating custom Detector Metadata

The scripts in the `python` folder use the tools in the detray python detector package to define various detector metadata headers in the `include` folder. A new generator file can be added here to the `python` folder, which will then be used during the cmake configuration step to generate the required detector metadata header for compilation of custom detector types (add `-DDETRAY_GENERATE_METADATA="my_metadata.py"` to the cmake setup)

To run any of the scripts manually, use the following:
```shell
source path-to-detray-build/python/setup.sh
python3 detray/detectors/python/odd_metadata.py
```

### Using generated metadata in client code

The generated headers are only available if they were requested in `DETRAY_GENERATE_METADATA`. Do not include them directly. Instead, include the aggregate header, which CMake always creates:

```cpp
#include <detray/detectors/generated_metadata.hpp>
```

For a generator `my_metadata.py`, which writes the header `my_metadata.hpp` with the struct `detray::my_metadata`, the aggregate header contains:

- `#define DETRAY_HAS_MY_METADATA`: Use this feature macro to guard code that needs this specific metadata type.
- `X(my)` in the X-macro `DETRAY_FOR_EACH_GENERATED_METADATA(X)`: Use it to write code for all generated metadata types without naming them, e.g.

```cpp
#define MAKE_DETECTOR_TYPE(NAME) \
  using NAME##_detector = detray::host::detector<detray::NAME##_metadata<algebra_t>>;
DETRAY_FOR_EACH_GENERATED_METADATA(MAKE_DETECTOR_TYPE)
#undef MAKE_DETECTOR_TYPE
```

Generated headers with the same name as a static header in `include/detray/detectors` are not included in the aggregate header.
