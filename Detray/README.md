# detray

Detray is a modern, C++20 header-only library providing a GPU-friendly tracking detector description using different linear algebra libraries. It follows the navigation and propagation concept of ACTS, but implementing a geometry using a flat memory layout and no abstract interfaces (virtual functions). A detray detector can therefore be constructed on the host and copied to an accelerator device in a straight-forward way.

With the geometry description comes a fully featured, GPU-ready track state propagation implementation in inhomogeneous magnetic fields (vector field description using [covfie](https://github.com/acts-project/covfie)), with track parameter covariance transport including material interactions.

## Requirements and Dependencies
#### OS & Compilers:

- The C++ compiler must support C++20
- The CUDA Toolkit version must be greater than major version 11

#### Dependencies:
- CMake (version >= 3.21)


## Getting started

The repository should build "out of the box", with standard
CMake build procedures.

```shell
git clone https://github.com/acts-project/acts.git
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -S Acts/Detray -B detray-build
cmake --build detray-build
```

For unit and integration tests using the *Open Data Detector* (ODD) solenoid field, a magnetic field map file in covfie format needs to be downloaded and the corresponding environment variable should be set to:
```shell
cd Detray/data
bash detray_data_get_files.sh
export DETRAY_BFIELD_FILE="${PWD}/odd-bfield_v0_9_0.cvf"
```

#### Build options

A number of cmake preset configurations are provided and can be listed by:
```shell
cmake -S Acts/Detray --list-presets
```
For a developer build, the `detray-dev-fp32` and `detray-dev-fp64` configurations are available (`fp`: floating point precision):
```shell
cmake -S Acts/Detray -B detray-build --preset detray-dev-fp32
```
The developer presets will build the components of detray that are most commonly used. The `prefetch` presets on the other hand will configure all possible dependencies, but not automatically trigger the build of the corresponding components. For example, in order to trigger the build of the unit tests with the `prefetch` preset, the corresponding option needs to be specified:
```shell
cmake -S Acts/Detray -B detray-build --preset detray-prefetch-fp32 \
-DDETRAY_ARRAY_PLUGIN=ON -DDETRAY_BUILD_UNITTESTS=ON
```
A full build, containing all components (e.g. tests and benchmarks), can be configured using the `detray-full-fp32` and `detray-full-fp64` presets.

The following cmake options are available and can also be specified explicitly for any preset:

| Option | Description | Default |
| --- | --- | --- |
| DETRAY_SET_LOGGING  | Set log level (NONE, WARN, INFO, VERBOSE, DEBUG) | INFO |
| DETRAY_BUILD_CUDA  | Build the CUDA sources included in detray | ON (if available) |
| DETRAY_BUILD_SYCL  | Build the SYCL sources included in detray | OFF |
| DETRAY_BUILD_TEST_UTILS  | Build the detray test utilities library (contains e.g. test detectors) | OFF |
| DETRAY_BUILD_UNITTESTS  | Build the detray unit tests | OFF |
| DETRAY_BUILD_INTEGRATIONTESTS  | Build the detray integration tests | OFF |
| DETRAY_BUILD_ALL_TESTS  | Build the detray unit and integration tests | OFF |
| DETRAY_BUILD_BENCHMARKS  | Build the detray benchmarks | OFF |
| DETRAY_BUILD_CLI_TOOLS  | Build the detray command line tools | OFF |
| DETRAY_BUILD_TUTORIALS  | Build the examples of detray | OFF |
| DETRAY_BUILD_PYTHON_BINDINGS  | Build python bindings for all enabled components | OFF |
| DETRAY_CUSTOM_SCALARTYPE | Floating point precision | float |
| DETRAY_GENERATE_METADATA | List of metadata generator scripts (separated by semicolon) | empty |
| DETRAY_EIGEN_PLUGIN | Build Eigen math plugin | OFF |
| DETRAY_FASTOR_PLUGIN | Build Fastor math plugin | OFF |
| DETRAY_SMATRIX_PLUGIN | Build ROOT/SMatrix math plugin | OFF |
| DETRAY_VC_AOS_PLUGIN | Build Vc based AoS math plugin | OFF |
| DETRAY_VC_SOA_PLUGIN | Build Vc based SoA math plugin (currently only supports the ray-surface intersectors) | OFF |
| DETRAY_SVG_DISPLAY | Build ActSVG display module | OFF |

## Tutorials

In the `tutorials` folder of the repository, there are a number of standalone executables that showcase different tasks:
- Building a detector (either a predefined one or designing a detector from scratch)
- Detector file IO
- Additional debugging options, like performing a ray scan or obtaining
extra navigation tracing information
- Moving a detector to device
- Host and device track propagation

In order to define a custom detector geometry type (called a detector 'metadata'), please follow the instructions in `Acts/Detray/detectors/README.md`.

Otherwise, the default detector metadata (`#include Acts/Detray/detectors/default_metadata.hpp`) can be used in most cases to define the detector type, however, incurring increased build times and likely also increased runtime of client algorithms.

## Detector Validation

Given a detray detector (and optionally also a grid and a material) json file, a number of validation test can be run from the command-line. For this, the library has to be built with the `-DDETRAY_BUILD_CLI_TOOLS=ON` option enabled. An example detector file can then be obtained using e.g.
```shell
detray-build/bin/detray_generate_toy_detector --write_material --write_grids
```
All of the validation tools presented in the following can also be run as part of a corresponding [python script](https://github.com/acts-project/Acts/Detray/tests/tools/python) which takes the same arguments and will automatically create plots from the collected data. However, this requires Python 3, pandas, SciPy and NumPy, as well as Matplotlib to be available.

The detector geometry can be visualized in SVG format with the following command:
```shell
detray-build/bin/detray_detector_display \
   --geometry_file  ./toy_detector/toy_detector_geometry.json
```
The tool can also display single volumes or surfaces, as well as the navigation grids and material maps (the corresponding json files need to loaded in this case). For an overview of all available options for the command-line tools add `--help`.

### Navigation Validation

In order to validate that the navigation works correctly in a given detector geometry, run the detector validation tool. It will first perform a consistency check on the detector, followed by a "ray scan" of the detector. The scan result will be compared to a full straight-line navigation run for every ray. After that, the navigation in a constant magnetic field of 2T is being tested in a similar fashion, using parameterized helix trajectories and a Newton-Raphson/Bisection algorithm to generate the truth intersections. For example:
```shell
detray-build/bin/detray_navigation_validation \
    --geometry_file ./toy_detector/toy_detector_geometry.json \
    --grid_file ./toy_detector/toy_detector_surface_grids.json \
    --search_window 3 3 --n_tracks 100 --pT_range 0.5 100
```
In case of failures, this command will give a detailed debug output in the form of a log file, as well as an SVG representation of the failed tracks. The grid file is optional, but will trigger the use of spatial grids as acceleration structures during the navigation run.

Note: The `search_window` option defines the size of lookup area of the grid acceleration structure and is therefore detector dependent! Use `--search_window 3 3` (or larger) for the *toy detector* and *wire chamber* example detectors and `--search_window 0 0` otherwise.

### Material Validation

This tool checks whether the navigator picks up the material correctly by comparing the material found during a ray scan with the material collected during navigation by a specialized actor:
```shell
detray-build/bin/detray_material_validation \
    --geometry_file ./toy_detector/toy_detector_geometry.json \
    --material_file ./toy_detector/toy_detector_homogeneous_material.json \
    --phi_steps 100 --eta_steps 100 --eta_range -4 4
```
Note: The correct material file must be loaded in addition to the geometry file!


## Benchmarks

A number of benchmarks exist, which are based on the google benchmark library, and can be run from command-line. For this, the `-DDETRAY_BUILD_BENCHMARKS=ON` and `-DDETRAY_BUILD_CLI_TOOLS=ON` flags need to be specified. Then pass the detray detector file(s) and additional options to the benchmark tools for the different hardware backends:
```shell
detray-build/bin/detray_propagation_benchmark_<backend>_<algebra> \
    --geometry_file ./toy_detector/toy_detector_geometry.json \
    --grid_file ./toy_detector/toy_detector_surface_grids.json \
    --material_file ./toy_detector/toy_detector_homogeneous_material.json \
    --sort_tracks --randomize_charge --eta_range -3 3 --pT_range 0.5 100 \
    --search_window 3 3 --covariance_transport --bknd_name [HW_BACKEND_NAME]
```
For every algebra-plugin that was built, a corresponding benchmark executable will be present. The CPU-backend benchmark is built by default and the CUDA-backend benchmark will be available if detray was built with CUDA enabled (`-DDETRAY_BUILD_CUDA=ON`). This executable can additionally be configured with any arguments targeted at [google benchmark](https://github.com/google/benchmark/blob/main/docs/user_guide.md).

If the data is dumped into json files using the options `--benchmark_repetitions=N` (N standing for the number of repetitions), `--benchmark_display_aggregates_only=true`, as well as `--benchmark_out_format=json` and `--benchmark_out=<some_file_name>.json`, it can afterwards be plotted with e.g.:
```shell
python3 detray/tests/tools/python/propagation_benchmarks.py \
    --geometry_file ./toy_detector/toy_detector_geometry.json \
    --data_files [FILES]...
```
