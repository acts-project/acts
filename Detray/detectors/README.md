## Generating custom Detector Metadata

The scripts in the `python` folder use the tools in the detray python detector package to define various detector metadata headers in the `include` folder. A new generator file can be added here to the `python` folder, which will then be used during the cmake configuration step to generate the required detector metadata header for compilation of custom detector types (add `-DDETRAY_GENERATE_METADATA="my_metadata.py"` to the cmake setup)

To run any of the scripts manually, use the following:
```shell
source path-to-detray-build/python/setup.sh
python3 detray/detectors/python/odd_metadata.py
```
